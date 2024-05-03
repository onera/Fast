        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)

        il    = l                     -  inck  
        ir    = l

        l200  = l  -  inci
        l202  = l  -  inci            -  inck  
        l102  = l  +  inci            -  inck  
        l100  = l  +  inci                    
        l022  = l             -  incj -  inck  
        l020  = l             -  incj  
        l012  = l             +  incj -  inck  
        l010  = l             +  incj  
CCCCCC
CCCCCC
CCCCCC   Facette K
CCCCCC
CCCCCC
      tiz  =  tcz
      tky  =  0.25*tcy
      tjx  =  0.25*tcx

      ! --- Gradients en U_k
      f1 = rop(il +v2)
      f2 = rop(ir +v2)
      fv = f1 + f2
      f3 = fv + rop(l202 +v2)+rop(l200 +v2)
      f4 = fv + rop(l102 +v2)+rop(l100 +v2)

      gradU_nx=(f4-f3)*tjx
      gradU_nz=(f2-f1)*tiz

      ! --- Gradients en V_k
      f1 = rop(il +v3)
      f2 = rop(ir +v3)
      fv = f1 + f2
      f5 = fv + rop(l022 +v3)+rop(l020 +v3)
      f6 = fv + rop(l012 +v3)+rop(l010 +v3)

      gradV_ny=(f6-f5)*tky
      gradV_nz=(f2-f1)*tiz

      ! --- Gradients en W_k
      f1 = rop(il +v4)
      f2 = rop(l  +v4)
      fv = f1 + f2
      f3 = fv + rop(l202 +v4)+rop(l200 +v4)
      f4 = fv + rop(l102 +v4)+rop(l100 +v4)
      f5 = fv + rop(l022 +v4)+rop(l020 +v4)
      f6 = fv + rop(l012 +v4)+rop(l010 +v4)

      gradW_nx=(f4-f3)*tjx
      gradW_ny=(f6-f5)*tky
      gradW_nz=(f2-f1)*tiz
        
      ! --- Gradients en T_k
      f1 =  rop(il +v5)
      f2 =  rop(ir +v5)

      gradT_nz=(f2-f1)*tiz

      !--- assemblage flux  K
      div = ( gradU_nx+gradV_ny+ gradW_nz)*cvisq

      f3 =  gradU_nz + gradW_nx
      f5 =  gradV_nz + gradW_ny
      f6 =  gradW_nz - div

#include  "FastS/Compute/mut_prandtltb_interface.for"

      fv     =  -   f3*tcz*xmutvol
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     =  -   f5*tcz*xmutvol
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     =  -2.*f6*tcz*xmutvol
      fv5    = fv5  + fv * (rop(ir+v4)+rop(il +v4))
      flu4   = flu4 + fv

      fv     = fv5*0.5  - gradT_nz*tcz *xktvol
      flu5   = flu5 + fv
