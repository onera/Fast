c***********************************************************************
c     $Date: 2010-01-28 16:22:02 +0100 (Thu, 28 Jan 2010) 
c     $ $Revision: 56 $ 
c     $Author: IvanMary $
c*****a*****************************************************************
      subroutine extrap_rk3localpara(param_int, ind_loop,rop) 
      
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

      INTEGER_E ind_loop(6), param_int(0:*)

      REAL_E rop(param_int(NDIMDX),param_int(NEQ))
 
C Var local
      INTEGER_E l,ijkm,im,jm,km,ldjr,i,j,k,ne,lij,neq,n,nistk,lstk
      INTEGER_E ni,nj,shift


#include "FastS/formule_param.h"

      ni = param_int(IJKV)     + 2*param_int(NIJK+3)
      nj = param_int(IJKV + 1) + 2*param_int(NIJK+3)

      neq = param_int(NEQ)
      
      !!! extrapolation dans la direction -i
      if (ind_loop(1).eq.1) then
          do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                  do i = -param_int(NIJK+3)+1, 0
                     
                     l  = inddm(i,j,k)
                     shift = 1 - i
                     rop(l,ne) = rop(l + shift , ne)

                   end do
                end do
              end do
            end do
       end if

      !!! extrapolation dans la direction i
      if (ind_loop(2).eq.param_int(IJKV)) then
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(3), ind_loop(4)
                  do i = ind_loop(2)+1, ind_loop(2) + param_int(NIJK+3)
                     
                     l  = inddm(i,j,k)
                  
                     
                     shift = ind_loop(2) - i
                     rop(l,ne) = rop(l + shift , ne)

                     
                   end do
                end do
              end do
            end do
        end if
       
                   
      !!! extrapolation dans la direction -j
      if (ind_loop(3).eq.1) then
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = -param_int(NIJK+3)+1, 0
                  do i = ind_loop(1),ind_loop(2)
                     


                     l  = inddm(i,j,k)
                     shift = 1 - j
                     rop(l,ne) = rop(l + shift*ni, ne)


                   end do
                end do
              end do
            end do
         end if

      !!! extrapolation dans la direction j
      if (ind_loop(4).eq.param_int(IJKV+1)) then
         do  ne=1,neq
             do  k = ind_loop(5), ind_loop(6)
               do  j = ind_loop(4)+1, ind_loop(4) + param_int(NIJK+3)
                  do i = ind_loop(1),ind_loop(2)
                     
                     l  = inddm(i,j,k)         
                     shift = ind_loop(4) - j 
                     rop(l,ne) = rop(l + shift*ni, ne)

                   end do
                end do
              end do
            end do
          end if

c$$$      !!! extrapolation dans la direction -k
c$$$      if (ind_loop(5).eq.1 .and. param_int(NIJK+3).ne.0) then
c$$$         do  ne=1,neq
c$$$             do  k = -param_int(NIJK+3)+1,0
c$$$               do  j = ind_loop(3), ind_loop(4)
c$$$                  do i = ind_loop(1),ind_loop(2)
c$$$                     
c$$$                     l  = inddm(i,j,k)
c$$$                     shift = 1 - k
c$$$                     rop(l,ne) = rop(l + shift*ni*nj, ne)
c$$$
c$$$                   end do
c$$$                end do
c$$$              end do
c$$$            end do
c$$$          end if
c$$$
c$$$
c$$$      !!! extrapolation dans la direction k
c$$$      if (ind_loop(5).eq.param_int(IJKV+2) .and. param_int(NIJK+3).ne.0) then
c$$$         do  ne=1,neq
c$$$             do  k = ind_loop(6)+1,ind_loop(6) + param_int(NIJK+3)
c$$$               do  j = ind_loop(3), ind_loop(4)
c$$$                  do i = ind_loop(1),ind_loop(2)
c$$$                     
c$$$                     l  = inddm(i,j,k)
c$$$                     shift = ind_loop(6) - k
c$$$                     rop(l,ne) = rop(l + shift*ni*nj, ne)
c$$$
c$$$                   end do
c$$$                end do
c$$$              end do
c$$$            end do
c$$$          end if             

      
      end
