c***********************************************************************
c     $Date: 2011-12-07 15:30:38 +0100 (mer 07 dec 2011) $
c     $Revision: 58 $
c     $Author: MarcTerracol $
c***********************************************************************
      subroutine synchro_omp_scater(param_int, ithread_io,
     &                          Nbre_thread_actif,
     &                          lok_shap, neq_lok,
     &                          ithread, topo_thread, thread_pos,
     &                          synchro_receive_th  , synchro_send_th,
     &                          ith  , jth  ,  kth  , ijkv_th,
     &                          size_cache,
     &                          ind_dm_thread,
     &                          lok, ctype)
c***********************************************************************
c_P                          O N E R A
c=======================================================================
      implicit none

      include 'omp_lib.h'

      character*7 ctype
      INTEGER_E ithread_io, neq_lok,
     & Nbre_thread_actif,
     & size_cache(3),
     & ithread, topo_thread(3), thread_pos(3),
     & synchro_send_th(3),synchro_receive_th(3),
     & ith, jth, kth, ijkv_th(3),
     & ind_dm_thread(6),
     & lok_shap(4),
     & lok( lok_shap(4), lok_shap(3), neq_lok, Nbre_thread_actif),
     & param_int(0:*)


#include "FastC/param_solver.h"

C Var loc  
      character*20 str
      logical ldir2,ldir3
      integer*8 compteur
      INTEGER_E l,lth, iverbs, ithread_rac,
     & ldir,lmod, no, ithread_i1,l0,
     & verrou_cachebloc(3), shift, overlap, i,
     & sens(3)

#define DEBUG 0

      if(ind_dm_thread(2).lt.ind_dm_thread(1)) return

      iverbs = 1

      ithread_i1 = ithread_io

      no = param_int (IO_THREAD )

      lth =1
      lmod =   mod(lth, neq_lok)
      if(lmod.eq.0) lmod = neq_lok

      ldir3 =.false.
      ldir2 =.false.
      if(ijkv_th(2).ne.1) ldir2 =.true.
      if(ijkv_th(3).ne.1) ldir3 =.true.


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

          Endif !Imin

          !!Imax
          If( ith.eq.ijkv_th(1)) then

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

          Endif!Jmin  

          
          !!Jmax
          If( jth.eq.ijkv_th(2) ) then

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
          Endif !Kmin

      ELSEIF(ctype(1:5).eq.'wait ') THEN  ! mise en attente du thread pour le calcul du sousdomaine ndsdm :


c      write(*,*)'topo', topo_thread
c      write(*,*)'pos', thread_pos
c      write(*,*)'ijk', ijkv_th
c      write(*,*)'synch', synchro_receive_th
          sens(1)=1
          sens(2)=1
          sens(3)=1
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
          Endif 

          !
          !
          !synchro Jmin
          !
          !
          If(jth.eq.1) Then

           !! Thread
           if(topo_thread(2).ne.1) then
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
           Endif

          !
          !
          !synchro Kmin: rien a faire
          !
          !


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

           Endif !Kmax

      ELSEIF(ctype(1:7).eq.'wait_lu') THEN  !mise en attente du thread pour LU ou mise a jour explithit

       write(*,*)'stop lu_wait'
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
