c***********************************************************************
c     $Date: 2011-12-07 15:30:38 +0100 (mer 07 dec 2011) $
c     $Revision: 58 $
c     $Author: MarcTerracol $
c***********************************************************************
      subroutine synchro_omp_scater(param_int, ithread_io,
     &                          lth,sens,lgo,lwait,Nbre_socket,
     &                          Nbre_thread_actif,thread_parsock, 
     &                          lok_sock_shap, lok_shap, neq_lok,
     &                          socket , topo_socket, socket_pos,
     &                          ithread, topo_thread, thread_pos,
     &                          synchro_receive_sock, synchro_send_sock,
     &                          synchro_receive_th  , synchro_send_th,
     &                          isock, jsock, ksock , ijkv_sock,
     &                          ith  , jth  ,  kth  , ijkv_th,
     &                          size_cache,
     &                          ind_dm_thread,
     &                          lok_new, lok_sock, lok, ctype)
c***********************************************************************
c_P                          O N E R A
c=======================================================================
      implicit none

      include 'omp_lib.h'

      character*7 ctype
      INTEGER_E ithread_io, lth, neq_lok,
     & Nbre_thread_actif,thread_parsock,lgo,lwait,Nbre_socket,
     & socket , topo_socket(3), socket_pos(3), sens(3), size_cache(3),
     & ithread, topo_thread(3), thread_pos(3),
     & synchro_send_sock(3),synchro_receive_sock(3),
     & synchro_send_th(3),synchro_receive_th(3),
     & isock, jsock,ksock, ith,jth,kth, ijkv_sock(3), ijkv_th(3),
     & ind_dm_thread(6), lok_new(Nbre_thread_actif),
     & lok_sock_shap(4),lok_shap(4),
     & lok_sock(lok_sock_shap(4), lok_sock_shap(3), Nbre_socket),
     & lok( lok_shap(4), lok_shap(3), neq_lok, Nbre_thread_actif),
     & param_int(0:*)


#include "FastC/param_solver.h"

C Var loc  
      character*20 str
      logical llu_wait, lwait_i, lwait_j, lwait_k,ldir2,ldir3
      integer*8 compteur
      INTEGER_E l, iverbs,ncell, blk_shift,ithread_rac,socket_rac,
     & ldir,lmod, no, imax,jmax,kmax,ithread_i1,itest,l0,
     & verrou_cachebloc(3),shift,overlap,i

#define DEBUG 0

      if(ind_dm_thread(2).lt.ind_dm_thread(1)) return

      iverbs = 1

      ithread_i1 = ithread_io

      no = param_int (IO_THREAD )
c      lth= isock+ (jsock-1)*ijkv_sock(1) 
c     &          + (ksock-1)*ijkv_sock(1)*ijkv_sock(2)

      lmod =   mod(lth, neq_lok)
      if(lmod.eq.0) lmod = neq_lok

      ldir3 =.false.
      ldir2 =.false.
      if(ijkv_th(2).ne.1) ldir2 =.true.
      if(ijkv_th(3).ne.1) ldir3 =.true.

      if(topo_thread(1).ne.1.and.ijkv_sock(1).ne.1) then
         lwait_i =.true.
         lwait_j =.false.
         lwait_k =.false.
      elseif(topo_thread(2).ne.1.and.ijkv_sock(2).ne.1) then
         lwait_i =.false.
         lwait_j =.true.
         lwait_k =.false.
      elseif(topo_thread(3).ne.1.and.ijkv_sock(3).ne.1) then 
         lwait_i =.false.
         lwait_j =.false.
         lwait_k =.true.
      else
         lwait_i =.false.
         lwait_j =.false.
         lwait_k =.false.
      endif


      IF(ctype(1:3).eq.'go ') THEN  !Feu vert pour sousdomaine ndsdm
       
          !!Imin
          If(ith.eq.1) then

             !!Thread
             if(synchro_send_th(1).eq.1
     &          .and.(thread_pos(1).ne.1.and.topo_thread(1).ne.1)) Then

               l   = kth + (jth-1)*ijkv_th(3) 
               ldir= 1

#include       "FastC/HPC_LAYER/flush.for"

               str='avgominI, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou(  lok(l, ldir, lmod, ithread) , 1 , compteur)

               str='  gominI, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"
            endif !thread        
            if(     isock.eq.1 .and. synchro_send_sock(1).eq.1
     &         .and.jth.eq.ijkv_th(2).and.kth.eq.ijkv_th(3)
     &         .and.thread_pos(2).eq.topo_thread(2)
     &         .and.thread_pos(3).eq.topo_thread(3)
     &         .and. thread_pos(1).eq.1                      ) Then 

               l   = ksock + (jsock-1)*ijkv_sock(3) 
               ldir= 1

#include       "FastC/HPC_LAYER/flush.for"

               str='avgominI, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou(  lok_sock(l, ldir, socket) , 1 , compteur)
               !call verrou(  lok_sock(l, ldir, ithread) , 1 , compteur)

               str='  gominI, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"
            endif !Sock
          Endif !Imin

          !!Imax
          If( ith.eq.ijkv_th(1)) then

            !thread go bloc
            if(lwait_i.and.jth*kth.eq.ijkv_th(1)*ijkv_th(2)
     &                .and.lgo.eq.1) then 

                l   = 1
                ldir= lok_shap(3)

#include        "FastC/HPC_LAYER/flush.for"

                 str='avgobloI, ith...kth='
#include         "FastC/HPC_LAYER/write_go.for"

                 call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)

                 str='  gobloI, ith...kth='
#include         "FastC/HPC_LAYER/write_go.for"
            endif

              !Thread coin (imax,Jmin)
            if (    topo_thread(1).ne.1 .and.jth.eq.1
     &          .and.thread_pos(1).ne.topo_thread(1) 
     &          .and.thread_pos(2).ne.1                )  Then 

                l   = kth 
                ldir= 2

#include        "FastC/HPC_LAYER/flush.for"

                str='avgomaxI, ith...kth='
#include        "FastC/HPC_LAYER/write_go.for"

                call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)
             
                str='  gomaxI, ith...kth='
#include        "FastC/HPC_LAYER/write_go.for"
            endif

            !Sock coin (imax,Jmin)
            if ( topo_socket(1).ne.1.and.jsock.eq.1.
     &                              .and.isock.eq.ijkv_sock(1)
     &                              .and.jth.eq.1
     &                              .and.kth.eq.ijkv_th(3)
     &                              .and.socket_pos(1).ne.topo_socket(1)
     &                              .and.thread_pos(1).eq.topo_thread(1)
     &                              .and.thread_pos(3).eq.topo_thread(3)
     &                              .and.thread_pos(2).eq.1     )  then 
                l   = ksock 
                ldir= 2
#include        "FastC/HPC_LAYER/flush.for"

                 str='avgomaxI, iso...kso='
#include         "FastC/HPC_LAYER/write_go.for"

                call verrou(  lok_sock(l, ldir, socket) , 1 , compteur)
                !call verrou(  lok_sock(l, ldir, ithread) , 1 , compteur)
                str='  gomaxI, iso...kso='
#include        "FastC/HPC_LAYER/write_go.for"
            endif
          Endif !Imax

          !!Jmin
          If(jth.eq.1) then

            !!Thread
            if( synchro_send_th(2).eq.1
     &         .and. (thread_pos(2).ne.1.and.topo_thread(2).ne.1) ) then

                l   = ith + (kth-1)*ijkv_th(1) 
                ldir= 1 + lok_shap(1)

#include        "FastC/HPC_LAYER/flush.for"

                str='avgominJ, ith...kth='
#include        "FastC/HPC_LAYER/write_go.for"

                call verrou( lok(l, ldir, lmod, ithread) , 1 , compteur)

                str='  gominJ, ith...kth='
#include        "FastC/HPC_LAYER/write_go.for"
            endif!tread

            if(   jsock.eq.1 .and. synchro_send_sock(2).eq.1
     &         .and.ith.eq.ijkv_th(1).and. kth.eq.ijkv_th(3)
     &          .and. thread_pos(1).eq.topo_thread(1)
     &          .and. thread_pos(3).eq.topo_thread(3)
     &         .and.thread_pos(2).eq.1)                      then

                l   = isock + (ksock-1)*ijkv_sock(1) 
                ldir= 1 + lok_sock_shap(1)
#include        "FastC/HPC_LAYER/flush.for"

                 str='avgominJ, iso...kso='
#include         "FastC/HPC_LAYER/write_go.for"

                call verrou(  lok_sock(l, ldir, socket) , 1 , compteur)
                !call verrou(  lok_sock(l, ldir, ithread) , 1 , compteur)
                str='  gominJ, iso...kso='
#include        "FastC/HPC_LAYER/write_go.for"
            endif!Socket
          Endif!Jmin  

          
          !!Jmax
          If( jth.eq.ijkv_th(2) ) then

            ! Thread go bloc
            if(lwait_j.and.ith*kth.eq.ijkv_th(1)*ijkv_th(3)
     &                .and.thread_pos(2).eq.topo_thread(2)
     &                .and.lgo.eq.1) then 

               l   = 1
               ldir= lok_shap(3)
#include       "FastC/HPC_LAYER/flush.for"

               str='avgobloJ, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)

               str='  gobloJ, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"
             endif

             !Thread coin (jmax,kmin)
             if( topo_thread(2).ne.1 .and.kth.eq.1
     &                            .and.thread_pos(2).ne.topo_thread(2) 
     &                            .and.thread_pos(3).ne.1 )  then  
               l   = ith 
               ldir= 2 + lok_shap(1)

#include       "FastC/HPC_LAYER/flush.for"

               str='avgomaxJ, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)

               str='  gomaxJ, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"
             endif!coin

             !Socket coin (jmax,kmin)
             if(     topo_socket(2).ne.1 .and.ksock.eq.1
     &          .and.jsock.eq.ijkv_sock(2)
     &          .and.ith.eq.ijkv_th(1).and.kth.eq.1
     &          .and.socket_pos(2).ne.topo_socket(2) 
     &          .and. thread_pos(1).eq.topo_thread(1)
     &          .and. thread_pos(2).eq.topo_thread(2)
     &          .and.thread_pos(3).eq.1 )                  then  

               l   = isock 
               ldir= 2 + lok_sock_shap(1)
#include       "FastC/HPC_LAYER/flush.for"

               str='avgomaxJ, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou(  lok_sock(l, ldir, socket) , 1 , compteur)
               !call verrou(  lok_sock(l, ldir, ithread) , 1 , compteur)
               str='  gomaxJ, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"
             endif!coin
          Endif !Jmax


          !!Kmin
          If(kth.eq.1) Then
             !!Thread
             If( synchro_send_th(3).eq.1
     &          .and. (thread_pos(3).ne.1.and.topo_thread(3).ne.1)) Then

               l   = ith + (jth-1)*ijkv_th(1) 
               ldir= 1 + lok_shap(1)+ lok_shap(2)

#include       "FastC/HPC_LAYER/flush.for"

                l0 = l+ (ldir-1)*lok_shap(4)
     &                + (lmod-1)* lok_shap(3)*lok_shap(4)
     &                + (ithread-1)*neq_lok* lok_shap(3)*lok_shap(4)
               str='avgominK, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)

               str='  gominK, ith...kth='
#include       "FastC/HPC_LAYER/write_go.for"
             endif
             !!Socket
             if(     ksock.eq.1 .and. synchro_send_sock(3).eq.1
     &          .and. ith.eq.ijkv_th(1).and.jth.eq.ijkv_th(2)
     &          .and. thread_pos(1).eq.topo_thread(1)
     &          .and. thread_pos(2).eq.topo_thread(2)
     &          .and. thread_pos(3).eq.1)                         Then

               l   = isock + (jsock-1)*ijkv_sock(1) 
               ldir= 1 + lok_sock_shap(1)+ lok_sock_shap(2)

#include       "FastC/HPC_LAYER/flush.for"

               str='avgominK, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"

               call verrou(  lok_sock(l, ldir, socket) , 1 , compteur)
               !call verrou(  lok_sock(l, ldir, ithread) , 1 , compteur)
               str='  gominK, iso...kso='
#include       "FastC/HPC_LAYER/write_go.for"
             endif
          Endif !Kmin

           !!Kmax Thread
           If(lwait_k.and.ith*jth.eq.ijkv_th(1)*ijkv_th(2)
     &               .and.kth.eq.ijkv_th(3)
     &               .and.thread_pos(3).eq.topo_thread(3)
     &               .and.lgo.eq.1) then 

             l   = 1
             ldir= lok_shap(3)

#include     "FastC/HPC_LAYER/flush.for"

             str='avgomaxK, ith...kth='
#include     "FastC/HPC_LAYER/write_go.for"

             call verrou( lok(l, ldir, lmod, ithread), 1 , compteur)

             str='  gomaxK, ith...kth='
#include     "FastC/HPC_LAYER/write_go.for"
           Endif


      if(lgo.eq.2.and.ith.eq.ijkv_th(1)
     &           .and.jth.eq.ijkv_th(2)
     &           .and.kth.eq.ijkv_th(3)) then

             l    =-1
             ldir =-1
             lmod =-1
             str='avlgonew, ith...kth='
#include     "FastC/HPC_LAYER/write_go.for"

             call verrou( lok_new(ithread), 1 , compteur)
      endif

      ELSEIF(ctype(1:5).eq.'wait ') THEN  ! mise en attente du thread pour le calcul du sousdomaine ndsdm :

          !Determination de la postion du verrou pour gerer
          !recouvrememnt ind_coe sur 4 maille
          ! ijkv_th   si size_cache>=overlap
          ! ijkv_th-1 si size_cache=3
          ! ijkv_th-2 si size_cache=2
          ! ijkv_th-4 si size_cache=1
          overlap = 4
          do i=1,3
            verrou_cachebloc(i)= ijkv_th(i)
            if (size_cache(i).lt.overlap) then
               shift = 1 + mod(overlap,size_cache(i))/size_cache(i)
               verrou_cachebloc(i)=  max( 1, verrou_cachebloc(i) -shift)
            endif
          enddo
 
          !
          !
          !synchro Imin
          !
          !
          !!Thread 
          If(topo_thread(1).ne.1.and.ith.eq.1) Then

            if(lwait_i.and.isock*jsock*ksock.ne.1
     &                .and.lwait  .eq.1
     &                .and.jth*kth.eq.1) then

              ithread_rac =  ithread -  1
              if(thread_pos(1).eq.1)ithread_rac=ithread+topo_thread(1)-1
          
               lmod =   mod(lth-1, neq_lok)
               if(lmod.eq.0) lmod = neq_lok
               ldir= lok_shap(3)
               l    = 1

               str='avwait Iblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"

               call verrou( lok(l,ldir,lmod, ithread_rac), 2, compteur)

               str='  wait Iblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
            endif


            !traitement du coin Imin, Jmax
            if(jth.eq.ijkv_th(2).and.synchro_receive_th(2).eq.1
     &                         .and.thread_pos(1).ne.1
     &                         .and.thread_pos(2).ne.topo_thread(2))then

               ithread_rac = ithread + topo_thread(1)*sens(2)-sens(1)
               l           = kth
               ldir        = 2
               lmod =   mod(lth, neq_lok)
               if(lmod.eq.0) lmod = neq_lok

               str='avt coinIblo:ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
               call verrou( lok(l,ldir, lmod, ithread_rac), 2, compteur)

               str='waitcoinIblo:ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
            endif
          Endif
          !!Socket coin Imin, Jmax
          itest= ith*kth*isock*synchro_receive_sock(2)*thread_pos(1)
          If(     jsock.eq.ijkv_sock(2).and.itest.eq.1
     &       .and.jth.eq.ijkv_th(2)
     &       .and.socket_pos(1) .ne.1 .and.topo_socket(1).ne.1
     &       .and.socket_pos(2) .ne.topo_socket(2)
     &       .and.thread_pos(1).eq.1
     &       .and.thread_pos(3).eq.1
     &       .and.thread_pos(2) .eq.topo_thread(2) ) then

               socket_rac = socket + topo_socket(1)
               ithread_rac= socket_rac
c               ithread_rac= ithread + (topo_socket(1)-1)*thread_parsock !shift socket
c     &                              - (topo_thread(2)-1)*topo_thread(1) !jpos =1 
c     &                              + (topo_thread(1)-1)                !ipos = topo_thread(1)
               l           = ksock
               ldir        = 2

               str='avt coinIblo:isokso='
#include       "FastC/HPC_LAYER/write_wait.for"
               call verrou(  lok_sock(l,ldir, socket_rac) , 1, compteur)
               !call verrou( lok_sock(l,ldir,  ithread_rac), 2, compteur)
               str='waitcoinIblo:isokso='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
          Endif


          !
          !
          !synchro Imax
          !
          !
          !Thread
          !if(ith.eq.ijkv_th(1).and.thread_pos(1).ne.topo_thread(1)
          If( ith.eq.verrou_cachebloc(1) ) Then

             !!Thread
             if(     thread_pos(1).ne.topo_thread(1)
     &          .and.synchro_receive_th(1).eq.1 ) Then

                lmod =   mod(lth, neq_lok)
                if(lmod.eq.0) lmod = neq_lok

                ithread_rac =  ithread +  sens(1)

                l    = kth + (jth-1)*ijkv_th(3) 

                ldir = 1

                str='avwait Imax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"

                call verrou(lok(l,ldir, lmod, ithread_rac), 2, compteur)

                str='  wait Imax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif
             itest = kth*jth*synchro_receive_sock(1)
            !Socket
             if(     isock.eq.ijkv_sock(1)
     &          .and.thread_pos(2).eq.1
     &          .and.thread_pos(3).eq.1
     &          .and.itest.eq.1.and.thread_pos(1).eq.topo_thread(1))then

                socket_rac =  socket +  1
                ithread_rac= socket_rac
c               ithread_rac=  ithread +  thread_parsock !shift socket
c     &                               - (topo_thread(1)-1) !ipos =1

                l    = ksock + (jsock-1)*ijkv_sock(3) 
                ldir = 1
                str='avwait Imax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"
                call verrou( lok_sock(l,ldir,  socket_rac), 2, compteur)
               !call verrou( lok_sock(l,ldir,  ithread_rac), 2, compteur)
                str='  wait Imax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif 
          Endif 

          !
          !
          !synchro Jmin
          !
          !
          If(jth.eq.1) Then

           !! Thread
           if(topo_thread(2).ne.1) then

            if(lwait_j .and.isock*jsock*ksock.ne.1
     &                 .and.ith*kth.eq.1
     &                 .and.lwait  .eq.1
     &                 .and.thread_pos(2).eq.1 ) then

               if(thread_pos(2).eq.1) then
                ithread_rac =ithread +topo_thread(1)*(topo_thread(2) -1)
               else
                ithread_rac =  ithread - topo_thread(1) 
               endif
               ldir = lok_shap(3)
               l    = 1
               lmod =   mod(lth-1, neq_lok)
               if(lmod.eq.0) lmod = neq_lok

               str='avwait Jblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"

               call verrou( lok(l,ldir,lmod, ithread_rac), 2, compteur)

               str='  wait Jblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
            endif
            !traitement du coin Jmin, Kmax
            if(kth.eq.ijkv_th(3).and.synchro_receive_th(3).eq.1
     &                         .and.thread_pos(2).ne.1
     &                         .and.thread_pos(3).ne.topo_thread(3))then

               ithread_rac = ithread + topo_thread(2)*sens(3)-sens(2)
               l           = ith
               ldir        = 2 + lok_shap(1) 
               lmod =   mod(lth, neq_lok)
               if(lmod.eq.0) lmod = neq_lok

               str='avt coinJblo:ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
               call verrou( lok(l,ldir, lmod, ithread_rac), 2, compteur)

               str='waitcoinJblo:ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
            endif
           Endif !topo_thread

           !!Socket coin Jmin, kmax
           itest= jth*ith*jsock*synchro_receive_sock(3)*thread_pos(2)
           if(   ksock.eq.ijkv_sock(3).and.itest.eq.1
     &       .and.kth.eq.ijkv_th(3)
     &       .and.socket_pos(2) .ne.1 .and.topo_socket(2).ne.1
     &       .and.socket_pos(3) .ne.topo_socket(3)
     &       .and.thread_pos(1).eq.1
     &       .and.thread_pos(2).eq.1
     &       .and.thread_pos(3) .eq.topo_thread(3) ) then

               socket_rac = socket + topo_socket(2)-1
               ithread_rac= socket_rac
c               ithread_rac= ithread + (topo_socket(2)-1)*thread_parsock !shift socket
c     &                              - (topo_thread(3)-1)*topo_thread(2) !kpos =1 
c     &                              + (topo_thread(2)-1)                !jpos = topo_thread(2)

               l           = isock
               ldir        = 2 + lok_sock_shap(1) 

               str='avt coinJblo:isokso='
#include       "FastC/HPC_LAYER/write_wait.for"
               call verrou( lok_sock(l,ldir, socket_rac),2,compteur)
               !call verrou( lok_sock(l,ldir, ithread_rac),2,compteur)
               str='waitcoinJblo:isokso='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
           endif
          Endif !Jmin



          !
          !
          !synchro Jmax
          !
          !
           !If( jth.eq.ijkv_th(2)) Then
           If( jth.eq.verrou_cachebloc(2) ) Then

             !!Thread
             if(      thread_pos(2).ne.topo_thread(2)
     &           .and.synchro_receive_th(2).eq.1      )  then

                ithread_rac = ithread + topo_thread(1)*sens(2)
                l           = ith + (kth-1)*ijkv_th(1) 

                ldir = 1 + lok_shap(1) 
                lmod =   mod(lth, neq_lok)
                if(lmod.eq.0) lmod = neq_lok

                str='avwait Jmax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"

                call verrou( lok(l, ldir, lmod, ithread_rac),2,compteur)

                str='  wait Jmax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif
        
             !!Socket
             itest = kth*ith*synchro_receive_sock(2)
             if(    jsock.eq.ijkv_sock(2) .and.itest.eq.1
     &          .and.thread_pos(1).eq.1
     &          .and.thread_pos(3).eq.1
     &          .and. thread_pos(2).eq.topo_thread(2) ) then

                socket_rac = socket + topo_socket(1)
                ithread_rac= socket_rac
c                ithread_rac = ithread +topo_socket(1)*thread_parsock
c     &                                -(topo_thread(2)-1)*topo_thread(1)!jpos =1
c
                l    = isock + (ksock-1)*ijkv_sock(1) 
                ldir = 1 + lok_sock_shap(1) 

                str='avwait Jmax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"

                call verrou( lok_sock(l, ldir, socket_rac),2,compteur)
                !call verrou( lok_sock(l, ldir, ithread_rac),2,compteur)

                str='  wait Jmax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif
           Endif

          !
          !
          !synchro Kmin
          !
          !
          If(lwait_k.and.isock*jsock*ksock.ne.1.and.ith*jth.eq.1
     &              .and.lwait  .eq.1
     &              .and. thread_pos(3).eq.1 
     &              .and. kth.eq.1) then

               if(thread_pos(3).eq.1) then
                ithread_rac =ithread +topo_thread(2)*(topo_thread(3)-1)
               else
                ithread_rac =  ithread - topo_thread(2) 
               endif

               ldir = lok_shap(3)
               l    = 1
               lmod =   mod(lth-1, neq_lok)
               if(lmod.eq.0) lmod = neq_lok

               str='avwait Kblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"

               call verrou( lok(l,ldir,lmod, ithread_rac), 2, compteur)

               str='  wait Kblo: ithkth='
#include       "FastC/HPC_LAYER/write_wait.for"
#include       "FastC/HPC_LAYER/flush.for"
           Endif

          !
          !
          !synchro Kmax
          !
          !
           !If( kth.eq.ijkv_th(3)) Then
           If( kth.eq.verrou_cachebloc(3) ) Then
             !!thread
             if(     thread_pos(3).ne.topo_thread(3)
     &          .and.synchro_receive_th(3).eq.1     )  then

                ithread_rac = ithread + topo_thread( 2 )*sens(3)
                l           = ith + (jth-1)*ijkv_th(1) 

                lmod =   mod(lth, neq_lok)
                if(lmod.eq.0) lmod = neq_lok
                ldir = 1 + lok_shap(1)+ lok_shap(2)

                l0 = l+ (ldir-1)*lok_shap(4)
     &                + (lmod-1)* lok_shap(3)*lok_shap(4)
     &                + (ithread_rac-1)*neq_lok* lok_shap(3)*lok_shap(4)
                str='avwait Kmax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"

                call verrou( lok(l, ldir, lmod, ithread_rac),2,compteur)

                str='  wait Kmax :ithkth='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif
             !!Socket
             itest = jth*ith*synchro_receive_sock(3)
             if(      ksock.eq.ijkv_sock(3).and. itest.eq.1
     &          .and.thread_pos(1).eq.1
     &          .and.thread_pos(2).eq.1
     &          .and. thread_pos(3).eq.topo_thread(3)           )then

                socket_rac = socket +topo_socket(2) 
                ithread_rac= socket_rac
c                ithread_rac = ithread +topo_socket(2)*thread_parsock
c     &                                -(topo_thread(3)-1)*topo_thread(2)!kpos =1

                l    = isock + (jsock-1)*ijkv_sock(1) 
                ldir = 1 + lok_sock_shap(1)+ lok_sock_shap(2)

                str='avwait Kmax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"
                call verrou( lok_sock(l, ldir, socket_rac),2,compteur)
                !call verrou( lok_sock(l, ldir, ithread_rac),2,compteur)
                str='  wait Kmax :isokso='
#include        "FastC/HPC_LAYER/write_wait.for"
#include        "FastC/HPC_LAYER/flush.for"
             endif
           Endif !Kmax

      if(lwait.eq.2.and.ith*jth*kth.eq.1
     &             .and.isock*jsock*ksock.ne.1) then

             l    =-1
             ldir =-1
             lmod =-1
             ithread_rac = ithread -1
             if(ithread_rac.lt.1+(socket-1)*thread_parsock) 
     &       ithread_rac =  socket*thread_parsock
             str='avlwait , ith...kth='
#include        "FastC/HPC_LAYER/write_wait.for"
      
             call verrou( lok_new(ithread_rac), 2 , compteur)
      endif

      ELSEIF(ctype(1:7).eq.'wait_lu') THEN  !mise en attente du thread pour LU ou mise a jour explithit

       write(*,*)'strop lu_wait'
       stop

C            compteur=0
C  13  continue
C            llu_wait = 0
C            do l=1,ndsdm_tot
C             call flush_integer( 1,lok(l))
C             if(lok(l).eq.0) llu_wait = 1
C            enddo     
C            compteur = compteur + 1
C!$          if(iverbs.ge.1.and.mod(compteur,100000000).eq.0)
C!$   &        write(*,*)'wait lu',ndsdm,compteur,ndsdm_tot   !ndom,compteur,nitcfg
C            if(llu_wait) goto 13
      ENDIF

      end
