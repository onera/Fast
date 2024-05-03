c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine conservrk3local4para(param_int,param_int2,coe,coe2,
     & param_real,ind_loop,ind_loop_,ind_loop__,rop,constk1,constk2,
     & tailleR,tailleD,nstep,dir,pt1,pt2,pt3,transfo,ratio1,ratio2,
     & ratio3)
    
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

      INTEGER_E idir,ind_loop(6),param_int(0:*),pos,tailleR,tailleD
      INTEGER_E param_int2(0:*), ind_loop_(6),dir,pt1,pt2,pt3,transfo(3)
      INTEGER_E ind_loop__(6),ratio1,ratio2,ratio3,nstep
 
      REAL_E param_real(0:*)
      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
      REAL_E constk1(tailleD,param_int(NEQ))
      REAL_E constk2(tailleR,param_int(NEQ))
      REAL_E coe(param_int(NDIMDX),param_int(NEQ_COE))
      REAL_E coe2(param_int2(NDIMDX),param_int2(NEQ_COE))
      REAL_E c1,c2
            
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk,cycl
      INTEGER_E nistk2,nistk3,nistkzv,nistk2zv,nistk3zv,lstkzv,lstkz
      INTEGER_E nistkz, nistk2z, nistk3z,posfluxzv,posflux,indrc
      INTEGER_E ind_loopz(6),ind_loopzv(6),i0,j0,k0,pt_pivot(3),indpivot
      REAL_E ro_old, u_old,v_old,w_old,t_old,roe_old,cv,cvinv
      INTEGER_E ijkdinf(3),ijkdinf_(3),nistk_,nistk2_,nistk3_,a,b,c,comp
      REAL_E bflu1,bflu2,bflu3,bflu4,bflu5
  
#include "FastS/formule_param.h"                          


      
      cycl = param_int(NSSITER)/param_int(LEVEL)
      neq=param_int(NEQ)

      cv = param_real(CVINF)
      cvinv=1.0/cv

    

      !print*, param_int(NDIMDX),param_int2(NDIMDX) 
     
       nistk  =(ind_loop_(2)- ind_loop_(1))+1
       nistk2 =(ind_loop_(4)- ind_loop_(3))+1
       nistk3 =(ind_loop_(6)- ind_loop_(5))+1

       nistk_  =(ind_loop__(2)- ind_loop__(1))+1
       nistk2_ =(ind_loop__(4)- ind_loop__(3))+1
       nistk3_ =(ind_loop__(6)- ind_loop__(5))+1


  
       indpivot =  (pt1 +1 - ind_loop__(1))
     &            +(pt2 - ind_loop__(3))*nistk_
     &            +(pt3 -ind_loop__(5))*nistk_*nistk2_


       !print*, 'indpivot= ',indpivot,pt1,pt2,pt3

       do i=1,3

          if (transfo(i)==1) then
             ijkdinf(i)=1
          else if (transfo(i)==-1) then
             ijkdinf(i)=-1
          else if (transfo(i)==2) then
             ijkdinf(i)= nistk_
          else if (transfo(i)==-2) then
             ijkdinf(i)= - nistk_
          else if (transfo(i)==3) then
             ijkdinf(i)= nistk_*nistk2_
          else if (transfo(i)==-3) then
             ijkdinf(i)= - nistk_*nistk2_
          end if

       end do
      
       !print*, 'indpivot= ',indpivot
       ijkdinf_ = ijkdinf 

       !print*, ijkdinf(1)
       !print*, ijkdinf(2)
       !print*, ijkdinf(3)

       ijkdinf(1)=ratio1*ijkdinf(1)
       ijkdinf(2)=ratio2*ijkdinf(2)
       ijkdinf(3)=ratio3*ijkdinf(3)

   
       do  k = ind_loop(5), ind_loop(6)
        do  j = ind_loop(3), ind_loop(4)
         do  i = ind_loop(1), ind_loop(2)                              
               
           l  = inddm(i,j,k)

           !if (i==ind_loop(1) .and. j== ind_loop(3)) then
           !   print*, 'l= ', l
           !end if


           lstk  =  (i+1 - ind_loop_(1))
     &              +(j - ind_loop_(3))*nistk
     &            +(k-ind_loop_(5))*nistk*nistk2

           bflu1=0.0
           bflu2=0.0
           bflu3=0.0
           bflu4=0.0
           bflu5=0.0
           
c$$$           comp = 0
c$$$
c$$$           do a=1,ratio1
c$$$              do b=1,ratio2
c$$$                 do c=1,ratio3
c$$$
c$$$                   indrc = indpivot
c$$$     &              + (i-ind_loop_(1))*ijkdinf(1)+float(a-1)*ijkdinf_(1)
c$$$     &              + (j-ind_loop_(3))*ijkdinf(2)+float(b-1)*ijkdinf_(2)
c$$$     &              + (k-ind_loop_(5))*ijkdinf(3)+float(c-1)*ijkdinf_(3)
c$$$
c$$$                    lstkzv = indrc
c$$$
c$$$                    bflu1 = bflu1 + constk2(lstkzv,1)
c$$$                    bflu2 = bflu2 + constk2(lstkzv,2)
c$$$                    bflu3 = bflu3 + constk2(lstkzv,3)
c$$$                    bflu4 = bflu4 + constk2(lstkzv,4)
c$$$                    bflu5 = bflu5 + constk2(lstkzv,5)
c$$$                    
c$$$                    comp = comp + 1
c$$$
c$$$                  !if (j == 1 .and. i.eq.1) then
c$$$                  !print*,"blu1,const2= ",bflu1, constk2(lstkzv,1) 
c$$$                  !end if
c$$$
c$$$                  end do
c$$$              end do
c$$$            end do





                   indrc = indpivot
     &              + (i-ind_loop_(1))*ijkdinf(1)!+float(a-1)*ijkdinf_(1)
     &              + (j-ind_loop_(3))*ijkdinf(2)!+float(b-1)*ijkdinf_(2)
     &              + (k-ind_loop_(5))*ijkdinf(3)!+float(c-1)*ijkdinf_(3)


          lstkzv = indrc

          !print*, l , lstk, lstkzv

          
          c1=coe(l,1)/dble(param_int(LEVEL))
          c2=coe(l,1)/dble(param_int2(LEVEL))
          !print*, 'param_int2(LEVEL)= ',param_int2(LEVEL) 

             ro_old = rop(l,1)
              u_old = rop(l,2)
              v_old = rop(l,3)
              w_old = rop(l,4)
              t_old = rop(l,5)
           
             bflu1 = constk2(lstkzv,1)
             rop(l,1) = rop(l,1)+c1*constk1(lstk,1)+c2*bflu1!c1*bflu1!constk2(lstkzv,1)

             !if (j == 1 .and. i.eq.1) then
             ! print*,"rop= ",rop(l,1),ro_old,c1*constk1(lstk,1),c2*bflu1 
             !end if

             !print*, c1*constk1(lstk,1), c2*bflu1
             
             bflu2 = constk2(lstkzv,2)
             rop(l,2) =(ro_old*rop(l,2) + 
     & c1*constk1(lstk,2) + c2*bflu2)/rop(l,1)

             bflu3 = constk2(lstkzv,3)
             rop(l,3) =(ro_old*rop(l,3) + 
     & c1*constk1(lstk,3) + c2*bflu3)/rop(l,1)

             bflu4 = constk2(lstkzv,4)
             rop(l,4) =(ro_old*rop(l,4) + 
     & c1*constk1(lstk,4) + c2*bflu4)/rop(l,1)


             bflu5 = constk2(lstkzv,5)

             roe_old = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))


              rop(l,5) = cvinv*( ( roe_old +
     & c1*constk1(lstk,5) + c2*bflu5)/                              
     & rop(l,1)  -  0.5*(rop(l,2)*rop(l,2)
     &                + rop(l,3)*rop(l,3)
     &                + rop(l,4)*rop(l,4)))   



       !print*,lstkzv
      !if (i==ind_loop(1) .and. j==200 ) then
      !   print*,c1*constk1(lstk,1),' ',c1*bflu1
      !end if

           end do
          end do
         end do
       



      end 
     
       
