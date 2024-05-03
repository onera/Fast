      !vitesse interface depuis etat G et D ordre1
      qp2 = rop(nm +v2) 
      qp3 = rop(nm +v3)
      qp4 = rop(nm +v4)
      qm2 = rop( l +v2) 
      qm3 = rop( l +v3)
      qm4 = rop( l +v4)

!In:  qpm(2-4) et qen      
!Out: qn1, qn2 (vitesse interface depuis etat G et D
#include "FastS/Compute/Vit_ent/qn_ale_3dhomo_i.for"

      qp1 = rop(nm +v1) 
      qm1 = rop( l +v1) 
      qm6 = rop(l +v6)*qm1*qn2
      qp6 = rop(nm+v6)*qp1*qn1

      !vitesse interface moyenne Roe depuis etat G et D ordre1
      qp1 = sqrt( qp1 )
      qm1 = sqrt( qm1 )
      r_1 = 1./(qp1 + qm1)
      qp2 = r_1*( qp1*qp2 + qm1*qm2)
      qp3 = r_1*( qp1*qp3 + qm1*qm3)
      qp4 = r_1*( qp1*qp4 + qm1*qm4)
!In:  qpm(2-4) et qen      
!Out: qn1, qn2 (vitesse interface depuis etat G et D
#include "FastS/Compute/Vit_ent/qn_ale_3dhomo_i.for"

      tdu  = max(c1/s_1, abs(qn1) )   ! s_1 = 1./Surf

      flu6 = 0.5*(qm6 + qp6 - tdu*( rop(l +v6)*rop(l +v1)
     &                             -rop(nm+v6)*rop(nm+v1)
     &                            )
     &          )

c             if(ndom.eq.1.or.ndom.eq.4) then
c               test1 =isnan(flu6)
c               if (test1) then
c                    write(*,*) 'flu6',l-lij+1,j,k
c                    stop
c               endif
c              endif

c      qm6  = qm*r2*qn2
c      qp6  = qp*r1*qn1

c      flu6 = 0.5*(qm6 + qp6 +sign(1.,flu1)*(qp6-qm6))
