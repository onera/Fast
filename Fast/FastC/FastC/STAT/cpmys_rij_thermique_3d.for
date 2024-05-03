          t          = rop(l,5)
          t_cn       = t*cn

          !  <rouT>,<rovT>,<rowT>
          moy(m, 14) = moy(m, 14) + rou_cn * t
          moy(m, 15) = moy(m, 15) + rov_cn * t
          moy(m, 16) = moy(m, 16) + row_cn * t
          !  <T>
          moy(m, 17) = moy(m, 17) + t_cn
          moy(m, 18) = moy(m, 18) + t_cn * t

          !   calcul de <nu>
          cmu        = c4/max(mu,1.e-30)
          nu_cn      = mu_cn*rho_1

            !dvardc(l,4,:)= Cp mu grad(T) / Pr  --> grad(T)
c          dvardc(l,4,1) = dvardc(l,4,1)*cmu
c          dvardc(l,4,2) = dvardc(l,4,2)*cmu
c          dvardc(l,4,3) = dvardc(l,4,3)*cmu

          moy(m, 1 + neq_tensrey) = moy(m, 1 + neq_tensrey) + nu_cn
          moy(m, 2 + neq_tensrey) = moy(m, 2 + neq_tensrey) 
c     &                          + ( dvardc(l,4,1)*dvardc(l,4,1)
c     &                             +dvardc(l,4,2)*dvardc(l,4,2) 
c     &                             +dvardc(l,4,3)*dvardc(l,4,3) )*nu_cn

