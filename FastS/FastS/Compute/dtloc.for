          !calcul coe(1) elsA like
          d      = 1./(tc2i*tc2i+tc2j*tc2j+tc2k*tc2k)

          !calcul coe(1) funk like
          !d      = 1./(max(tc2i,tc2j,tc2k))
          !d      = d*d

          dtvis  = 0.5*detj*d*gam1_1/xmut(l)*r

          sp     = sqrt(u*u+v*v+w*w) + c
          dtconv = sqrt(d)/sp
 
          dt     = param_real(CFL)*min(dtconv,dtvis)

          coe(l,1)=dt

