       nd_rdm = 1

       do while( (no_lu(nd_rdm).ne.no_rdm.and.
     &            nd_rdm.lt.nisdom_residu(nitcfg)) )
          nd_rdm = nd_rdm + 1
       enddo
       
       if(no_lu(nd_rdm).eq.no_rdm) then
         lcomput=.true.
       else
         lcomput=.false.
       endif

       IF(lcomput) THEN

          xinterm = 0.

          !no_rdm: mono indice parcours sousbloc
          k_lu=1+(no_rdm-1)/(ijkv_lu(1)*ijkv_lu(2))
          j_lu=1+(no_rdm-1-ijkv_lu(1)*ijkv_lu(2)*(k_lu-1))/ijkv_lu(1)
          i_lu=no_rdm-(j_lu-1)*ijkv_lu(1)-(k_lu-1)*ijkv_lu(1)*ijkv_lu(2)

          ind_loop_lu(1) = 1 + (i_lu-1)*size_ssdom(1)
          ind_loop_lu(3) = 1 + (j_lu-1)*size_ssdom(2)
          ind_loop_lu(5) = 1 + (k_lu-1)*size_ssdom(3)

          ind_loop_lu(2) = i_lu*size_ssdom(1)
          if(i_lu.eq.ijkv_lu(1)) ind_loop_lu(2) = ijkv(1)
          ind_loop_lu(4) = j_lu*size_ssdom(2)
          if(j_lu.eq.ijkv_lu(2)) ind_loop_lu(4) = ijkv(2)
          ind_loop_lu(6) = k_lu*size_ssdom(3)
          if(k_lu.eq.ijkv_lu(3)) ind_loop_lu(6) = ijkv(3)

          if (iflw.eq.3.and.les.eq.0) then


#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
            DO l=1,ijkv(1)*ijkv(2)*ijkv(3)

              k=1+(l-1)/(ijkv(1)*ijkv(2))
              j=1+(l-1-ijkv(1)*ijkv(2)*(k-1))/ijkv(1)
              i=l-(j-1)*ijkv(1)-(k-1)*ijkv(1)*ijkv(2)
#else
            do k = ind_loop_lu(5), ind_loop_lu(6)
            do j = ind_loop_lu(3), ind_loop_lu(4)
            do i = ind_loop_lu(1), ind_loop_lu(2)   
#endif
              l = inddm(i,j,k)
   
              xp=abs(drodm(l,1))
              if (xp.gt.xinterm(1,2)) xinterm(1,2)=xp
              xinterm(1,1) = xinterm(1,1)+drodm(l,1)*drodm(l,1)
        
              xp=abs(drodm(l,2))
              if(xp.gt.xinterm(2,2)) xinterm(2,2)=xp
              xinterm(2,1) = xinterm(2,1)+drodm(l,2)*drodm(l,2)
   
              xp=abs(drodm(l,3))
              if(xp.gt.xinterm(3,2))  xinterm(3,2)=xp
              xinterm(3,1) = xinterm(3,1)+drodm(l,3)*drodm(l,3)
   
              xp=abs(drodm(l,4))
              if(xp.gt.xinterm(4,2))  xinterm(4,2)=xp
              xinterm(4,1) = xinterm(4,1)+drodm(l,4)*drodm(l,4)

              xp=abs(drodm(l,5))
              if(xp.gt.xinterm(5,2)) xinterm(5,2)=xp
              xinterm(5,1) = xinterm(5,1)+drodm(l,5)*drodm(l,5)

              xp=abs(drodm(l,6))
              if(xp.gt.xinterm(6,2)) xinterm(6,2)=xp
              xinterm(6,1) = xinterm(6,1)+drodm(l,6)*drodm(l,6)
#ifndef E_SCALAR_COMPUTER
            enddo
#else
            enddo
            enddo
            enddo
#endif

          else

#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
            DO l=1,ijkv(1)*ijkv(2)*ijkv(3)

              k=1+(l-1)/(ijkv(1)*ijkv(2))
              j=1+(l-1-ijkv(1)*ijkv(2)*(k-1))/ijkv(1)
              i=l-(j-1)*ijkv(1)-(k-1)*ijkv(1)*ijkv(2)
#else
            do k = ind_loop_lu(5), ind_loop_lu(6)
            do j = ind_loop_lu(3), ind_loop_lu(4)
            do i = ind_loop_lu(1), ind_loop_lu(2)   
#endif
              l = inddm(i,j,k)

              xp=abs(drodm(l,1))
              if (xp.gt.xinterm(1,2)) xinterm(1,2)=xp
C              if (xp.gt.xinterm(1,2)) then
C                xinterm(1,2)=xp
C                ii = i
C                jj = j
C                kk = k
C              endif
              xinterm(1,1) = xinterm(1,1)+drodm(l,1)*drodm(l,1)

              xp=abs(drodm(l,2))
              if(xp.gt.xinterm(2,2)) xinterm(2,2)=xp
              xinterm(2,1) = xinterm(2,1)+drodm(l,2)*drodm(l,2)

              xp=abs(drodm(l,3))
              if(xp.gt.xinterm(3,2))  xinterm(3,2)=xp
              xinterm(3,1) = xinterm(3,1)+drodm(l,3)*drodm(l,3)

              xp=abs(drodm(l,4))
              if(xp.gt.xinterm(4,2))  xinterm(4,2)=xp
              xinterm(4,1) = xinterm(4,1)+drodm(l,4)*drodm(l,4)

              xp=abs(drodm(l,5))
              if(xp.gt.xinterm(5,2)) xinterm(5,2)=xp
              xinterm(5,1) = xinterm(5,1)+drodm(l,5)*drodm(l,5)
#ifndef E_SCALAR_COMPUTER
            enddo
#else
            enddo
            enddo
            enddo
#endif

          endif!test si SA/DES ou autre
         
          it_lu_ssdom(no_rdm) = it_lu_ssdom(no_rdm) + 1
          do ne=1,neq
            rdm_sdm(it_bloc, ne, 1, no_rdm) = xinterm(ne,1)
            rdm_sdm(it_bloc, ne, 2, no_rdm) = xinterm(ne,2)
          enddo

          xro1  = LOG10(max(rdm_sdm(it_bloc,1,2,no_rdm),cut0x))
     &          - LOG10(max(rdm_sdm(1      ,1,2,no_rdm),cut0x))
          xrou1 = LOG10(max(rdm_sdm(it_bloc,2,2,no_rdm),cut0x))
     &          - LOG10(max(rdm_sdm(1      ,2,2,no_rdm),cut0x)) 
          xrov1 = LOG10(max(rdm_sdm(it_bloc,3,2,no_rdm),cut0x))
     &          - LOG10(max(rdm_sdm(1      ,3,2,no_rdm),cut0x)) 
          xrow1 = LOG10(max(rdm_sdm(it_bloc,4,2,no_rdm),cut0x)) 
     &          - LOG10(max(rdm_sdm(1      ,4,2,no_rdm),cut0x))
          xroe1 = LOG10(max(rdm_sdm(it_bloc,5,2,no_rdm),cut0x))
     &          - LOG10(max(rdm_sdm(1      ,5,2,no_rdm),cut0x))

          if(nijk(5).ne.0) then
            rmax = (xrou1+xrov1+xrow1+xro1+xroe1)*0.2
          else
            rmax = (xrou1+xrov1+xro1+xroe1)*0.25
          endif

C          if(neq.eq.6) then
C          xroe1 = LOG10(max(rdm_sdm(it_bloc,6,2,no_rdm),cut0x))
C     &          - LOG10(max(rdm_sdm(1 ,6,2,no_rdm),cut0x))
C          rmax = (5.*rmax+ xroe1)/6.
C          endif

          if(rmax.gt.conv_loo) then

            it_temp_ssdom(no_rdm) = it_lu_ssdom(no_rdm)


            !Si convergence ricrac, on ajoute une iter de gras
            if(it_temp_ssdom(no_rdm).eq.it_target_ssdom(no_rdm)
     &         .and.nssiter.eq.nitcfg) then

               it_temp_ssdom(no_rdm)=it_temp_ssdom(no_rdm)+1
            endif
          endif

        ELSE! sous_domaine pas calcule a cette iteration

          it = it_bloc
          do ne=1,neq
            rdm_sdm(it,ne,1,no_rdm)= rdm_sdm(it-1,ne,1,no_rdm)
            rdm_sdm(it,ne,2,no_rdm)= rdm_sdm(it-1,ne,2,no_rdm)
          enddo

        ENDIF

       if(nitcfg.eq.nssiter) then

        it_target_old(no_rdm)  = it_target_ssdom(no_rdm)
        it_target_ssdom(no_rdm)=max(3, it_temp_ssdom(no_rdm)+1)

c      if(ndom.eq.0.and.ithread.eq.1)
c     & write(*,*)'tgt old',it_target_ssdom(no_rdm),
c     & it_target_old(no_rdm),no_rdm
       endif
