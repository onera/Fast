c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine conservrk3local2(idir,param_int,param_int2,coe,coe2,
     & param_real,ind_loop,rop,constk1,constk2,taille,nstep)
    
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    extrapolation ordere zero cell fictive
c
c     VAL
c_V    Optimisation NEC
c
c     COM
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir,ind_loop(6),param_int(0:*),pos,taille,nstep
      INTEGER_E param_int2(0:*)
 
      REAL_E param_real(0:*)
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E constk1(taille/param_int(NEQ),param_int(NEQ))
      REAL_E constk2(taille/param_int(NEQ),param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E coe2(param_int2(NDIMDX),param_int2(NEQ_COE))
      REAL_E c1,c2
            
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,cycl
      INTEGER_E nistk2,nistk3,nistkzv,nistk2zv,nistk3zv,lstkzv,lstkz
      INTEGER_E nistkz, nistk2z, nistk3z,posfluxzv,posflux
      INTEGER_E ind_loopz(6),ind_loopzv(6)
      REAL_E ro_old, u_old,v_old,w_old,t_old,roe_old,cv,cvinv

  
#include "FastS/formule_param.h"                          


      
      cycl = param_int(NSSITER)/param_int(LEVEL)
      neq=param_int(NEQ)

      cv = param_real(CVINF)
      cvinv=1.0/cv


    

      !print*, param_int(NDIMDX),param_int2(NDIMDX) 
     
       nistk = (ind_loop(2)- ind_loop(1))+1
       nistk2 =(ind_loop(4)- ind_loop(3))+1
       nistk3 =(ind_loop(6)- ind_loop(5))+1
  

       ind_loopzv(3) = 1
       ind_loopzv(4) = param_int2(IJKV+1)
       ind_loopzv(5) = 1
       ind_loopzv(6) = param_int2(IJKV+2)

       if (ind_loop(2) == 1) then
          ind_loopzv(1)=param_int2(IJKV)+1
          ind_loopzv(2)=param_int2(IJKV)+1
          posflux = 1
          posfluxzv = 0
       else 
          ind_loopzv(1)=1!param_int(IJKV)+1
          ind_loopzv(2)=1!param_int(IJKV)+1
          posflux=0
          posfluxzv=1
       end if

       !print*,'posflux= ', posflux,'posfluxzv= ', posfluxzv

       nistkzv  =(ind_loopzv(2)- ind_loopzv(1))+1
       nistk2zv =(ind_loopzv(4)- ind_loopzv(3))+1
       nistk3zv =(ind_loopzv(6)- ind_loopzv(5))+1



      
       do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)                              
               
           l  = inddm(i,j,k)


           lstk  =  (i+1 - ind_loop(1))
     &              +(j - ind_loop(3))*nistk
     &            +(k-ind_loop(5))*nistk*nistk2
     &         + posflux*nistk*nistk2*nistk3



          ! lstkzv  =  (i+1 - ind_loopzv(1))
     &    !          +(j - ind_loopzv(3))*nistkzv
     &    !        +(k-ind_loopzv(5))*nistkzv*nistk2zv
     &    !     + posfluxzv*nistkzv*nistk2zv*nistk3zv

          lstkzv  =  (i+1 - ind_loop(1))
     &              +(j - ind_loop(3))*nistk
     &            +(k-ind_loop(5))*nistk*nistk2
     &         + posfluxzv*nistk*nistk2*nistk3

          
          c1=coe(l,1)/dble(param_int(LEVEL))
          c2=coe(l,1)/dble(param_int2(LEVEL))

             ro_old = rop(l,1)
              u_old = rop(l,2)
              v_old = rop(l,3)
              w_old = rop(l,4)
              t_old = rop(l,5)
           

             rop(l,1) = rop(l,1)+c1*constk1(lstk,1)+c2*constk2(lstkzv,1)

             rop(l,2) =(ro_old*rop(l,2) + 
     & c1*constk1(lstk,2) + c2*constk2(lstkzv,2))/rop(l,1)

             rop(l,3) =(ro_old*rop(l,3) + 
     & c1*constk1(lstk,3) + c2*constk2(lstkzv,3))/rop(l,1)

             rop(l,4) =(ro_old*rop(l,4) + 
     & c1*constk1(lstk,4) + c2*constk2(lstkzv,4))/rop(l,1)




             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))


              rop(l,5) = cvinv*( ( roe_old +
     & c1*constk1(lstk,5) + c2*constk2(lstkzv,5))/                              
     & rop(l,1)  -  0.5*(rop(l,2)*rop(l,2)
     &                + rop(l,3)*rop(l,3)
     &                + rop(l,4)*rop(l,4)))   



       !print*,lstkzv
      !if (j==85.and.ind_loop(1)==param_int(IJKV).and.cycl==16) then
      !   print*,rop(l,1),constk1(lstk,1),constk2(lstkzv,1),c1,c2
      !end if

           end do
          end do
         end do
       



      end 
     
       
