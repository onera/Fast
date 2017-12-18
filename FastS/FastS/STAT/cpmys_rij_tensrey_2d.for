          cmu        = 1./max(mu,1.e-30)

c          du1dx1     = dvardc(l,1,1)*cmu
c          du1dx2     = dvardc(l,1,2)*cmu
c          du2dx1     = dvardc(l,2,1)*cmu
c          du2dx2     = dvardc(l,2,2)*cmu

          dukdxk     = du1dx1+du2dx2

          s11        = du1dx1+du1dx1
          s22        = du2dx2+du2dx2
          s12        = du1dx2+du2dx1

c            <u>,<v>,<w>
          moy(m, eq0 + 1) = moy(m, eq0 + 1) + u1 * cn
          moy(m, eq0 + 2) = moy(m, eq0 + 2) + u2 * cn

c             terme II --------------------------------------------
          moy(m, eq0 + 4) = moy(m, eq0 + 4) + u1* u1 *rou_cn
          moy(m, eq0 + 5) = moy(m, eq0 + 5) + u2* u2 *rov_cn
          moy(m, eq0 + 7) = moy(m, eq0 + 7) + u1* u2 *rou_cn
          moy(m, eq0 + 9) = moy(m, eq0 + 9) + u2* u2 *rou_cn
c             
c             terme IIIa --------------------------------------------
          moy(m, eq0 +14) = moy(m, eq0 +14) + mu_cn   *dukdxk
          moy(m, eq0 +15) = moy(m, eq0 +15) + mu_cn*u1*dukdxk
          moy(m, eq0 +16) = moy(m, eq0 +16) + mu_cn*u2*dukdxk

c             terme IIIb --------------------------------------------
          moy(m, eq0 +18) = moy(m, eq0 +18) + mu_cn*( u1*s11 + u1*s11 ) ! 111
          moy(m, eq0 +19) = moy(m, eq0 +19) + mu_cn*( u2*s22 + u2*s22 ) ! 222
          moy(m, eq0 +21) = moy(m, eq0 +21) + mu_cn*( u1*s12 + u1*s12 ) ! 112
          moy(m, eq0 +23) = moy(m, eq0 +23) + mu_cn*( u2*s12 + u2*s12 ) ! 221
          moy(m, eq0 +28) = moy(m, eq0 +28) + mu_cn*( u2*s12 + u1*s22 ) ! 122
          moy(m, eq0 +29) = moy(m, eq0 +29) + mu_cn*( u2*s11 + u1*s12 ) ! 121
                                
          moy(m, eq0 +36) = moy(m, eq0 +36) + mu_cn*s11 ! 11
          moy(m, eq0 +37) = moy(m, eq0 +37) + mu_cn*s22 ! 22
          moy(m, eq0 +39) = moy(m, eq0 +39) + mu_cn*s12 ! 12

c             terme IV --------------------------------------------
          moy(m, eq0 +42) = moy(m, eq0 +42) + p_cn*(u1)
          moy(m, eq0 +43) = moy(m, eq0 +43) + p_cn*(u2)

c             terme VI ------------------------------------------- 
          moy(m, eq0 +45) = moy(m, eq0 +45) + p_cn*s11 ! 11
          moy(m, eq0 +46) = moy(m, eq0 +46) + p_cn*s22 ! 22
          moy(m, eq0 +48) = moy(m, eq0 +48) + p_cn*s12 ! 12
          moy(m, eq0 +51) = moy(m, eq0 +51) + cn  *s11 ! 11
          moy(m, eq0 +52) = moy(m, eq0 +52) + cn  *s22 ! 22
          moy(m, eq0 +54) = moy(m, eq0 +54) + cn  *s12 ! 12
                                 
c             terme VIIa --------------------------------------------
          moy(m, eq0 +57) = moy(m, eq0 +57) + mu_cn*s11*dukdxk ! 11
          moy(m, eq0 +58) = moy(m, eq0 +58) + mu_cn*s22*dukdxk ! 22
          moy(m, eq0 +60) = moy(m, eq0 +60) + mu_cn*s12*dukdxk ! 12
                                 
c             terme VIIb --------------------------------------------
          moy(m, eq0 +63)= moy(m, eq0 +63)+mu_cn*(du1dx1*s11+du1dx1*s11)! 111
          moy(m, eq0 +64)= moy(m, eq0 +64)+mu_cn*(du2dx2*s22+du2dx2*s22)! 222
          moy(m, eq0 +66)= moy(m, eq0 +66)+mu_cn*(du1dx2*s12+du1dx2*s12)! 112
          moy(m, eq0 +68)= moy(m, eq0 +68)+mu_cn*(du2dx1*s12+du2dx1*s12)! 221
          moy(m, eq0 +73)= moy(m, eq0 +73)+mu_cn*(du2dx2*s12+du1dx2*s22)! 122
          moy(m, eq0 +74)= moy(m, eq0 +74)+mu_cn*(du2dx1*s11+du1dx1*s12)! 121
