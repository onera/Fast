c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine indice_boucle_scater(ndo, ibloc,jbloc,kbloc,ithread,no,
     &                                topo_s, size_thread, ijkv_thread,
     &                                ijkvloc, thread_pos,ind_dm_socket,
     &                                ind_dm_omp )

      implicit none

      INTEGER_E ndo, ibloc,jbloc,kbloc,
     & topo_s(3), size_thread(3), ijkv_thread(3), ijkvloc(3),
     & thread_pos(3),ind_dm_socket(6),
     & ind_dm_omp(6)  ,ithread,no

C Var loc 
      INTEGER_E i,c4,Jmin, itest,jdeb,jv_remain,size_loc

      do i =1,3

         c4=0
         if(topo_s(i).ne.1) c4=1

           Jmin = size_thread(i)*topo_s(i)
           if(i.eq.1) then 
            itest = ibloc
           elseif(i.eq.2) then
            itest = jbloc
           else
            itest = kbloc
           endif

           if(itest.ne.ijkv_thread(i).or.mod(ijkvloc(i),Jmin).eq.0) then

               ind_dm_omp(2*i-1)=ind_dm_socket(2*i-1)+(itest -1)*Jmin
     &                          +(thread_pos(i)-1)*size_thread(i)*c4

               ind_dm_omp(2*i) =ind_dm_omp(2*i-1) + size_thread(i)-1

c            if(ithread.eq.no.and.jbloc.eq.1) then
c             write(*,*)'dom0',ind_dm_omp(2*i-1),ind_dm_omp(2*i),i,kbloc
c            endif

           else   
               jdeb = ind_dm_socket(2*i-1) +(itest-1)*Jmin

               jv_remain = ind_dm_socket(2*i)- jdeb + 1
               size_loc  = jv_remain/topo_s(i)


c            if(ithread.eq.no.and.jbloc.eq.1) then
c             write(*,*)'indi',jdeb,jv_remain,size_loc,i,kbloc
c            endif
             
            if(thread_pos(i).le.mod(jv_remain,topo_s(i))) then

                ind_dm_omp(2*i-1)=jdeb+(thread_pos(i)-1)*(size_loc+1)*c4
                ind_dm_omp(2*i  )= ind_dm_omp(2*i-1) + size_loc

c            if(ithread.eq.no.and.jbloc.eq.1) then
c             write(*,*)'dom1',ind_dm_omp(2*i-1),ind_dm_omp(2*i),i,kbloc
c            endif

            else
                ind_dm_omp(2*i-1)= jdeb +(thread_pos(i)-1)*size_loc*c4
     &                                  + mod(jv_remain,topo_s(i))
                ind_dm_omp(2*i)  = ind_dm_omp(2*i-1) + size_loc-1

c            if(ithread.eq.no.and.jbloc.eq.1) then
c             write(*,*)'dom2',ind_dm_omp(2*i-1),ind_dm_omp(2*i),i,kbloc
c            endif

            endif
         endif
      enddo  ! direction i,j,k

      end
