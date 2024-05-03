          !calcul coe(1) elsA like
          d      = 1./(tc2i*tc2i+tc2j*tc2j+tc2k*tc2k)

          dtvis  = 0.5*detj*d*gam1_1/xmut(l)*r

          sp     = sqrt(u*u+v*v+w*w) + c
          dtconv = sqrt(d)/sp
 
          dt     = param_real(CFL)*min(dtconv,dtvis)

          coe(l,1)=dt*MIN(cellN(l),2.-cellN(l))
