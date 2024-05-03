          cmu        = 1./max(mu,1.e-30)

c          du1dx1     = dvardc(l,1,1)*cmu
c          du1dx2     = dvardc(l,1,2)*cmu
c          du1dx3     = dvardc(l,1,3)*cmu
c          du2dx1     = dvardc(l,2,1)*cmu
c          du2dx2     = dvardc(l,2,2)*cmu
c          du2dx3     = dvardc(l,2,3)*cmu
c          du3dx1     = dvardc(l,3,1)*cmu
c          du3dx2     = dvardc(l,3,2)*cmu
c          du3dx3     = dvardc(l,3,3)*cmu

          dukdxk     = du1dx1+du2dx2+du3dx3

          s11        = du1dx1+du1dx1
          s22        = du2dx2+du2dx2
          s33        = du3dx3+du3dx3
          s12        = du1dx2+du2dx1
          s13        = du1dx3+du3dx1
          s23        = du2dx3+du3dx2

c            <u>,<v>,<w>
          moy(m, eq0 + 1) = moy(m, eq0 + 1) + u1 * cn
          moy(m, eq0 + 2) = moy(m, eq0 + 2) + u2 * cn
          moy(m, eq0 + 3) = moy(m, eq0 + 3) + u3 * cn

c             terme II --------------------------------------------
          moy(m, eq0 + 4) = moy(m, eq0 + 4) + u1* u1 *rou_cn
          moy(m, eq0 + 5) = moy(m, eq0 + 5) + u2* u2 *rov_cn
          moy(m, eq0 + 6) = moy(m, eq0 + 6) + u3* u3 *row_cn
          moy(m, eq0 + 7) = moy(m, eq0 + 7) + u1* u2 *rou_cn
          moy(m, eq0 + 8) = moy(m, eq0 + 8) + u1* u3 *rou_cn
          moy(m, eq0 + 9) = moy(m, eq0 + 9) + u2* u2 *rou_cn
          moy(m, eq0 +10) = moy(m, eq0 +10) + u2* u2 *row_cn
          moy(m, eq0 +11) = moy(m, eq0 +11) + u3* u3 *rou_cn
          moy(m, eq0 +12) = moy(m, eq0 +12) + u3* u3 *rov_cn
          moy(m, eq0 +13) = moy(m, eq0 +13) + u1* u2 *row_cn
c             
c             terme IIIa --------------------------------------------
          moy(m, eq0 +14) = moy(m, eq0 +14) + mu_cn   *dukdxk
          moy(m, eq0 +15) = moy(m, eq0 +15) + mu_cn*u1*dukdxk
          moy(m, eq0 +16) = moy(m, eq0 +16) + mu_cn*u2*dukdxk
          moy(m, eq0 +17) = moy(m, eq0 +17) + mu_cn*u3*dukdxk

c             terme IIIb --------------------------------------------
          moy(m, eq0 +18) = moy(m, eq0 +18) + mu_cn*( u1*s11 + u1*s11 ) ! 111
          moy(m, eq0 +19) = moy(m, eq0 +19) + mu_cn*( u2*s22 + u2*s22 ) ! 222
          moy(m, eq0 +20) = moy(m, eq0 +20) + mu_cn*( u3*s33 + u3*s33 ) ! 333
          moy(m, eq0 +21) = moy(m, eq0 +21) + mu_cn*( u1*s12 + u1*s12 ) ! 112
          moy(m, eq0 +22) = moy(m, eq0 +22) + mu_cn*( u1*s13 + u1*s13 ) ! 113
          moy(m, eq0 +23) = moy(m, eq0 +23) + mu_cn*( u2*s12 + u2*s12 ) ! 221
          moy(m, eq0 +24) = moy(m, eq0 +24) + mu_cn*( u2*s23 + u2*s23 ) ! 223
          moy(m, eq0 +25) = moy(m, eq0 +25) + mu_cn*( u3*s13 + u3*s13 ) ! 331
          moy(m, eq0 +26) = moy(m, eq0 +26) + mu_cn*( u3*s23 + u3*s23 ) ! 332
          moy(m, eq0 +27) = moy(m, eq0 +27) + mu_cn*( u2*s13 + u1*s23 ) ! 123
          moy(m, eq0 +28) = moy(m, eq0 +28) + mu_cn*( u2*s12 + u1*s22 ) ! 122
          moy(m, eq0 +29) = moy(m, eq0 +29) + mu_cn*( u2*s11 + u1*s12 ) ! 121
          moy(m, eq0 +30) = moy(m, eq0 +30) + mu_cn*( u3*s13 + u1*s33 ) ! 133
          moy(m, eq0 +31) = moy(m, eq0 +31) + mu_cn*( u3*s12 + u1*s23 ) ! 132
          moy(m, eq0 +32) = moy(m, eq0 +32) + mu_cn*( u3*s11 + u1*s13 ) ! 131
          moy(m, eq0 +33) = moy(m, eq0 +33) + mu_cn*( u3*s12 + u2*s13 ) ! 231
          moy(m, eq0 +34) = moy(m, eq0 +34) + mu_cn*( u3*s22 + u2*s23 ) ! 232
          moy(m, eq0 +35) = moy(m, eq0 +35) + mu_cn*( u3*s23 + u2*s33 ) ! 233
                                
          moy(m, eq0 +36) = moy(m, eq0 +36) + mu_cn*s11 ! 11
          moy(m, eq0 +37) = moy(m, eq0 +37) + mu_cn*s22 ! 22
          moy(m, eq0 +38) = moy(m, eq0 +38) + mu_cn*s33 ! 33
          moy(m, eq0 +39) = moy(m, eq0 +39) + mu_cn*s12 ! 12
          moy(m, eq0 +40) = moy(m, eq0 +40) + mu_cn*s13 ! 13
          moy(m, eq0 +41) = moy(m, eq0 +41) + mu_cn*s23 ! 23

c             terme IV --------------------------------------------
          moy(m, eq0 +42) = moy(m, eq0 +42) + p_cn*(u1)
          moy(m, eq0 +43) = moy(m, eq0 +43) + p_cn*(u2)
          moy(m, eq0 +44) = moy(m, eq0 +44) + p_cn*(u3)

c             terme VI ------------------------------------------- 
          moy(m, eq0 +45) = moy(m, eq0 +45) + p_cn*s11 ! 11
          moy(m, eq0 +46) = moy(m, eq0 +46) + p_cn*s22 ! 22
          moy(m, eq0 +47) = moy(m, eq0 +47) + p_cn*s33 ! 33
          moy(m, eq0 +48) = moy(m, eq0 +48) + p_cn*s12 ! 12
          moy(m, eq0 +49) = moy(m, eq0 +49) + p_cn*s13 ! 13
          moy(m, eq0 +50) = moy(m, eq0 +50) + p_cn*s23 ! 23
          moy(m, eq0 +51) = moy(m, eq0 +51) + cn  *s11 ! 11
          moy(m, eq0 +52) = moy(m, eq0 +52) + cn  *s22 ! 22
          moy(m, eq0 +53) = moy(m, eq0 +53) + cn  *s33 ! 33
          moy(m, eq0 +54) = moy(m, eq0 +54) + cn  *s12 ! 12
          moy(m, eq0 +55) = moy(m, eq0 +55) + cn  *s13 ! 13
          moy(m, eq0 +56) = moy(m, eq0 +56) + cn  *s23 ! 23
                                 
c             terme VIIa --------------------------------------------
          moy(m, eq0 +57) = moy(m, eq0 +57) + mu_cn*s11*dukdxk ! 11
          moy(m, eq0 +58) = moy(m, eq0 +58) + mu_cn*s22*dukdxk ! 22
          moy(m, eq0 +59) = moy(m, eq0 +59) + mu_cn*s33*dukdxk ! 33
          moy(m, eq0 +60) = moy(m, eq0 +60) + mu_cn*s12*dukdxk ! 12
          moy(m, eq0 +61) = moy(m, eq0 +61) + mu_cn*s13*dukdxk ! 13
          moy(m, eq0 +62) = moy(m, eq0 +62) + mu_cn*s23*dukdxk ! 23
                                 
c             terme VIIb --------------------------------------------
          moy(m, eq0 +63)= moy(m, eq0 +63)+mu_cn*(du1dx1*s11+du1dx1*s11)! 111
          moy(m, eq0 +64)= moy(m, eq0 +64)+mu_cn*(du2dx2*s22+du2dx2*s22)! 222
          moy(m, eq0 +65)= moy(m, eq0 +65)+mu_cn*(du3dx3*s33+du3dx3*s33)! 333
          moy(m, eq0 +66)= moy(m, eq0 +66)+mu_cn*(du1dx2*s12+du1dx2*s12)! 112
          moy(m, eq0 +67)= moy(m, eq0 +67)+mu_cn*(du1dx3*s13+du1dx3*s13)! 113
          moy(m, eq0 +68)= moy(m, eq0 +68)+mu_cn*(du2dx1*s12+du2dx1*s12)! 221
          moy(m, eq0 +69)= moy(m, eq0 +69)+mu_cn*(du2dx3*s23+du2dx3*s23)! 223
          moy(m, eq0 +70)= moy(m, eq0 +70)+mu_cn*(du3dx1*s13+du3dx1*s13)! 331
          moy(m, eq0 +71)= moy(m, eq0 +71)+mu_cn*(du3dx2*s23+du3dx2*s23)! 332
          moy(m, eq0 +72)= moy(m, eq0 +72)+mu_cn*(du2dx3*s13+du1dx3*s23)! 123
          moy(m, eq0 +73)= moy(m, eq0 +73)+mu_cn*(du2dx2*s12+du1dx2*s22)! 122
          moy(m, eq0 +74)= moy(m, eq0 +74)+mu_cn*(du2dx1*s11+du1dx1*s12)! 121
          moy(m, eq0 +75)= moy(m, eq0 +75)+mu_cn*(du3dx3*s13+du1dx3*s33)! 133
          moy(m, eq0 +76)= moy(m, eq0 +76)+mu_cn*(du3dx2*s12+du1dx2*s23)! 132
          moy(m, eq0 +77)= moy(m, eq0 +77)+mu_cn*(du3dx1*s11+du1dx1*s13)! 131
          moy(m, eq0 +78)= moy(m, eq0 +78)+mu_cn*(du3dx1*s12+du2dx1*s13)! 231
          moy(m, eq0 +79)= moy(m, eq0 +79)+mu_cn*(du3dx2*s22+du2dx2*s23)! 232
          moy(m, eq0 +80)= moy(m, eq0 +80)+mu_cn*(du3dx3*s23+du2dx3*s33)! 233
