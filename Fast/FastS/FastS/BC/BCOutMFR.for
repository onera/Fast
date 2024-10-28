CC    IN : rop
CC    OUT: roi,ui,vi,wi,Ti,pi (etat intern)
CC    OUT: ro0,u0,v0,w0,T0,pi (etat extern)
CC    OUT: tcx,tcy,tcz, tnx,tny,tnz
#include  "FastS/BC/linearization_state.for"
CC    In : tnx,tny,tnz, c_ale, ventijk
CC    OUT: qen (vitesse entrainement)
#include  "FastS/BC/qn_ale.for"
C ...  Implementation elsA
      roc0 = SQRT(ro0*gam*p0)

      c0_2 = roc0/ro0
      c0_2 = c0_2*c0_2
c ... Solving ro value using C0 and C+ characteristic curves
c     Solving ax**2+bx+c = 0
      coefa = c0_2
      coefb = - roc0*(u0*tnx + v0*tny + w0*tnz - qen) - c0_2*roi
      coefc = roc0*qp(li)
      delta = coefb*coefb-4*coefa*coefc

      rop(l,1) = (-coefb+SQRT(delta))/(2*coefa)
      dp = c0_2*(rop(l,1) - ro0)
      dqn = dp/roc0
      rop(l,2) = ui - dqn*tnx
      rop(l,3) = vi - dqn*tny
      rop(l,4) = wi - dqn*tnz
      rop(l,5) = (dp+pi)/(rop(l,1)*rgp)
