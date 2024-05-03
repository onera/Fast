c***********************************************************************
c     $Date: 2011-10-10 16:10:53 +0200 (lun 10 oct 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cptaylor(ndom, neq_grad,
     &                    param_int, ind_loop,
     &                    rop  , dvardc, 
     &                    tke  , enst , compteur,vol )
c***********************************************************************
c_U   USER : C. laurent
c
c     ACT
c_A    calcul grandeur moyenne pour bilan rij
c_A    commande de calcul des statistiques 
c
c     VAL
c
c     INP
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom, neq_grad, ind_loop(6), param_int(0:*)
c
      REAL_E    rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E vol( param_int(NDIMDX_MTR) )
      REAL_E dvardc( param_int(NDIMDX) * neq_grad*3 )
      REAL_E tke, enst,compteur

c Var loc
      INTEGER_E i,j,k,l,lij, v1,v2,v3,v4,v5,ltij,lt,lvo,
     & Vdudx,Vdudy,Vdudz, Vdvdx,Vdvdy,Vdvdz,Vdwdx,Vdwdy,Vdwdz

      REAL_E    u1,u2,u3,w1,w2,w3,cmu,rho,vref

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)

      Vdudx = 0
      Vdudy =   param_int(NDIMDX)
      Vdudz = 2*param_int(NDIMDX)
      Vdvdx = 3*param_int(NDIMDX)
      Vdvdy = 4*param_int(NDIMDX)
      Vdvdz = 5*param_int(NDIMDX)
      Vdwdx = 6*param_int(NDIMDX)
      Vdwdy = 7*param_int(NDIMDX)
      Vdwdz = 8*param_int(NDIMDX)




c$$$      compteur = compteur +  (ind_loop(2) -ind_loop(1)+1)*
c$$$     &                       (ind_loop(4) -ind_loop(3)+1)*
c$$$     &                       (ind_loop(6) -ind_loop(5)+1)

      if(param_int(ITYPZONE).eq.2) then ! zone cartesienne

        DO k = ind_loop(5), ind_loop(6)
        DO j = ind_loop(3), ind_loop(4)
#include "FastS/POST/loopI_begin.for"

          compteur = compteur + vol(1)

          !print*, compteur , vol(lt)

          rho        = rop(l+ V1)
          u1         = rop(l+ V2)
          u2         = rop(l+ V3) 
          u3         = rop(l+ V4)

          tke        = tke + 0.5*rho*(u1*u1+u2*u2+u3*u3)*vol(1)

          w1         = (dvardc(l + Vdwdy) - dvardc(l + Vdvdz))
          w2         = (dvardc(l + Vdudz) - dvardc(l + Vdwdx))
          w3         = (dvardc(l + Vdvdx) - dvardc(l + Vdudy))

          enst       =enst + 0.5*rho*(w1*w1+w2*w2+w3*w3)*vol(1)

        end do
        end do
        end do

      else

        DO k = ind_loop(5), ind_loop(6)
        DO j = ind_loop(3), ind_loop(4)
#include "FastS/POST/loopI_begin.for"

          compteur = compteur + vol(lt)

          rho        = rop(l+ V1)
          u1         = rop(l+ V2)
          u2         = rop(l+ V3) 
          u3         = rop(l+ V4)

          tke        = tke + 0.5*rho*(u1*u1+u2*u2+u3*u3)*vol(lt)

          w1         = (dvardc(l + Vdwdy) - dvardc(l + Vdvdz))
          w2         = (dvardc(l + Vdudz) - dvardc(l + Vdwdx))
          w3         = (dvardc(l + Vdvdx) - dvardc(l + Vdudy))

          enst       =enst + 0.5*rho*(w1*w1+w2*w2+w3*w3)*vol(lt)

        end do
        end do
        end do
      endif

      end
