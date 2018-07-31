***********************************************************************
c     $Date: 2018-11-04 13:25:50 +0100 (Thu, 04 Nov 2018) $
c     $Revision: 64 $
c     $Author: Thibaut $
c***********************************************************************
      subroutine dp_dw_vect(param_int, param_real, ind_loop, rop,
     &                      vectin, vectout, size)
c***********************************************************************
c_U   USER : GUEGAN
c
c     ACT
c_A    Jacobienne dp_dw: 
c      p               : (density,Velocity,Temperature, nutilde)
c      w               : (density,Momentum,TotalEnergy, Densitynutilde)
c
c     VAL
c
c     I/O
c_/    vectin    : vecteur kryloc IN
c_/    vectout   : vecteur kryloc OUT
c
c***********************************************************************      
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E param_int(0:*), ind_loop(6), size
      REAL_E param_real(0:*)
      REAL_E     rop( param_int(NDIMDX) , param_int(NEQ) )
      REAL_E  vectin( size              , param_int(NEQ) )
      REAL_E vectout( param_int(NDIMDX) , param_int(NEQ) )

      INTEGER_E i, k, j, l, i_2,j_2,k_2,v1,v2,v3,v4,ls,
     & v5,i_size,j_size,indlu
      REAL_E cvinv, roinv, ke, 
     &     dP2_dW1, dP2_dW2, dP3_dW1, dP3_dW3, dP4_dW1, dP4_dW4, 
     &     dP5_dW1, dP5_dW2, dP5_dW3, dP5_dW4, dP5_dW5,
     &     dP6_dW1, dP6_dW6

#include "FastS/formule_param.h"

      indlu(i_2,j_2,k_2) = 1 + (i_2+v1) + (j_2+v2)*v4 + (k_2+v3)*v5

      cvinv = 1./param_real(CVINF)

      if (param_int(NB_RELAX) .GE. 2) then

         i_size = ind_loop(2) - ind_loop(1) + 1 +
     &        2 * param_int(NIJK + 3) !taille de la fenetre + ghostcells
         j_size = ind_loop(4) - ind_loop(3) + 1 +
     &        2 * param_int(NIJK + 3)

         v1 = param_int(NIJK+3) - ind_loop(1)
         v2 = param_int(NIJK+3) - ind_loop(3)
         v3 = param_int(NIJK+4) - ind_loop(5)
         v4 = i_size
         v5 = i_size * j_size

      else

         v1 = param_int(NIJK+3) - 1
         v2 = param_int(NIJK+3) - 1
         v3 = param_int(NIJK+4) - 1
         v4 = param_int(NIJK)
         v5 = param_int(NIJK) * param_int(NIJK+1)

      endif

      !zone 2D
      IF(param_int(ITYPZONE).eq.3) THEN

        if(param_int(NEQ).eq.5) then

           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, 1)
              ls = indlu(i, j, 1)

               ke   = 0.5 * (rop(l,2)*rop(l,2) +rop(l,3)*rop(l,3))

               roinv   = 1./rop(l,1)

               dP2_dW1 = - rop(l,2)*roinv
               dP2_dW2 = roinv

               dP3_dW1 = - rop(l,3)*roinv
               dP3_dW3 = roinv

               dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
               dP5_dW2 = dP2_dW1*cvinv
               dP5_dW3 = dP3_dW1*cvinv
               dP5_dW5 = roinv*cvinv 

               vectout(l,1) = vectin(ls,1)
               vectout(l,2) = dP2_dW1*vectin(ls,1)+ dP2_dW2*vectin(ls,2)
               vectout(l,3) = dP3_dW1*vectin(ls,1)+ dP3_dW3*vectin(ls,3)
               vectout(l,4) = 0.
               vectout(l,5) = dP5_dW1*vectin(ls,1)+ dP5_dW2*vectin(ls,2)
     &                       +dP5_dW3*vectin(ls,3)+ dP5_dW5*vectin(ls,5)

           enddo
           enddo

        else
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, 1)
              ls = indlu(i, j, 1)

              ke   = 0.5 * ( rop(l,2)*rop(l,2) +rop(l,3)*rop(l,3))

              roinv   = 1./rop(l,1)

              dP2_dW2 = roinv
              dP2_dW1 = - rop(l,2)*roinv

              dP3_dW1 = - rop(l,3)*roinv
              dP3_dW3 = roinv

              dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
              dP5_dW2 = dP2_dW1*cvinv
              dP5_dW3 = dP3_dW1*cvinv
              dP5_dW5 = roinv*cvinv 

              dP6_dW1 = - rop(l,6)*roinv
              dP6_dW6 = roinv

              vectout(l,1) = vectin(ls,1)
              vectout(l,2) = dP2_dW1 *vectin(ls,1)+ dP2_dW2*vectin(ls,2)
              vectout(l,3) = dP3_dW1 *vectin(ls,1)+ dP3_dW3*vectin(ls,3)
              vectout(l,4) = 0.
              vectout(l,5) = dP5_dW1 *vectin(ls,1)+ dP5_dW2*vectin(ls,2)
     &             +dP5_dW3 *vectin(ls,3) + dP5_dW5*vectin(ls,5)
              vectout(l,6) = dP6_dW1 *vectin(ls,1)+ dP6_dW6*vectin(ls,6)

           enddo
           enddo
        endif

      !!
      !!zone 3d
      !!
      ELSE


        if(param_int(NEQ).eq.5) then
           do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, k)
              ls = indlu(i, j, k)

              ke   = 0.5*(  rop(l,2)*rop(l,2)
     &             + rop(l,3)*rop(l,3)
     &             + rop(l,4)*rop(l,4))

              roinv   = 1./rop(l,1)

              dP2_dW1 = - rop(l,2)*roinv
              dP2_dW2 = roinv

              dP3_dW1 = - rop(l,3)*roinv
              dP3_dW3 = roinv

              dP4_dW1 = - rop(l,4)*roinv
              dP4_dW4 = roinv

              dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
              dP5_dW2 = dP2_dW1*cvinv
              dP5_dW3 = dP3_dW1*cvinv
              dP5_dW4 = dP4_dW1*cvinv
              dP5_dW5 = roinv*cvinv 

              vectout(l,1) = vectin(ls,1)
              vectout(l,2) = dP2_dW1*vectin(ls,1) + dP2_dW2*vectin(ls,2)
              vectout(l,3) = dP3_dW1*vectin(ls,1) + dP3_dW3*vectin(ls,3)
              vectout(l,4) = dP3_dW1*vectin(ls,1) + dP4_dW4*vectin(ls,4)
              vectout(l,5) = dP5_dW1*vectin(ls,1) + dP5_dW2*vectin(ls,2)
     &             +dP5_dW3*vectin(ls,3) + dP5_dW4*vectin(ls,4)
     &             +dP5_dW5*vectin(ls,5)
            enddo
            enddo
            enddo
        else
           do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

              l = inddm(i, j, k)
              ls = indlu(i, j, k)

              ke   = 0.5*(  rop(l,2)*rop(l,2)
     &             + rop(l,3)*rop(l,3)
     &             + rop(l,4)*rop(l,4))

              roinv   = 1./rop(l,1)

              dP2_dW1 = - rop(l,2)*roinv
              dP2_dW2 = roinv

              dP3_dW1 = - rop(l,3)*roinv
              dP3_dW3 = roinv

              dP4_dW1 = - rop(l,4)*roinv
              dP4_dW4 = roinv

              dP5_dW1 = roinv*( ke*Cvinv - rop(l, 5) )
              dP5_dW2 = dP2_dW1*cvinv
              dP5_dW3 = dP3_dW1*cvinv
              dP5_dW4 = dP4_dW1*cvinv
              dP5_dW5 = roinv*cvinv 

              vectout(l,1) = vectin(ls,1)
              vectout(l,2) = dP2_dW1*vectin(ls,1) + dP2_dW2*vectin(ls,2)
              vectout(l,3) = dP3_dW1*vectin(ls,1) + dP3_dW3*vectin(ls,3)
              vectout(l,4) = dP3_dW1*vectin(ls,1) + dP4_dW4*vectin(ls,4)
              vectout(l,5) = dP5_dW1*vectin(ls,1) + dP5_dW2*vectin(ls,2)
     &             +dP5_dW3*vectin(ls,3) + dP5_dW4*vectin(ls,4)
     &             +dP5_dW5*vectin(ls,5)
              vectout(l,6) = dP6_dW1*vectin(ls,1) + dP6_dW6*vectin(ls,6)

            enddo
            enddo
            enddo
        endif

      ENDIF

      end
