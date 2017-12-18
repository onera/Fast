c***********************************************************************
c     $Date: 2011-12-07 15:30:38 +0100 (mer 07 déc 2011) $
c     $Revision: 56 $
c     $Author: MarcTerracol $
c***********************************************************************
      subroutine crsdm_scater(ndo,ith,jth,kth, topo_th,
     &                 size_cache, synchro_receive, synchro_send,
     &                 ind_dm_glob, ind_dm, ijkv_sdm )
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Creation de la partition d'un doamine en sous-domaine pour // openmp
c
c     INP    ndo
c
c     OUT    common paralle.h: ijkdeb_sdm
c
c=======================================================================
      implicit none

      include "omp_lib.h"

      INTEGER_E ndo, ith,jth,kth, topo_th(3),
     & size_cache(3), synchro_receive(3),synchro_send(3),
     & ind_dm_glob(6), ind_dm(6), ijkv_sdm(3)

C Var loc
      INTEGER_E l,i,j,k,ivloc,iverbs

      iverbs =0
 
        ! Faut il decouper dans la direction IJK pour gerer synchro des threads
        do i=1,3

         synchro_receive(i) = 0
         synchro_send   (i) = 0
         if( ind_dm_glob(2*i  ).ne. ind_dm(2*i).and.topo_th(i).ne.1) 
     &          synchro_receive(i)= 1
         if( ind_dm_glob(2*i-1).ne. ind_dm(2*i-1).and.topo_th(i).ne.1) 
     &       synchro_send(i)   = 1

         ivloc       = ind_dm(2*i)- ind_dm(2*i-1) + 1
        
         ijkv_sdm(i) = max(1, ivloc / size_cache(i) )

       enddo
   
C!$         if(iverbs.ge.4.and.OMP_get_thread_num().eq.0) 
C!$   &      write(*,'(8i5)')ndo, ijkdeb_sdm(1:6,ndsdm ),ndsdm


C!$         if(iverbs.ge.4.and.OMP_get_thread_num().eq.0)
C!$   &       write(*,*)'ndo et IJjkv_sdm(3)',ndo, ijkv_sdm

      end
