        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        il    = l             -  incj  
        ir    = l

        l022  = l             -  incj   -  inck  
        l002  = l                       -  inck  
        l021  = l             -  incj   +  inck  
        l001  = l                       +  inck  
        l220  = l  -  inci    -  incj  
        l200  = l  -  inci
        l120  = l  +  inci    -  incj  
        l100  = l  +  inci                    
CCCCCC
CCCCCC
CCCCCC   Facette J
CCCCCC
CCCCCC
      !Calcul des vecteurs surfaces J
      tiy  = tcy
      tjz  = .25 * tcz
      tkx  = .25 * tcx

      ! --- Gradients en U_j (i,j,k), (i,j-1,k), (i,j-1,k-1), (i,j,k-1)
      f1 =rop(il +v2)
      f2 =rop(ir +v2)
      fv = f1 + f2
      f5 = fv + rop(l220 +v2)+rop(l200 +v2)
      f6 = fv + rop(l120 +v2)+rop(l100 +v2)

      gradU_nx=(f6-f5)*tkx
      gradU_ny=(f2-f1)*tiy

      ! --- Gradients en V_j
      f1 =rop(il +v3)
      f2 =rop(ir +v3)
      fv = f1 + f2
      f3 = fv + rop(l022 +v3)+rop(l002 +v3)
      f4 = fv + rop(l021 +v3)+rop(l001 +v3)
      f5 = fv + rop(l220 +v3)+rop(l200 +v3)
      f6 = fv + rop(l120 +v3)+rop(l100 +v3)

      gradV_nx=(f6-f5)*tkx
      gradV_ny=(f2-f1)*tiy
      gradV_nz=(f4-f3)*tjz
      ! --- Gradients en W_j
      f1 =rop(il +v4)
      f2 =rop(ir +v4)
      fv = f1 + f2
      f3 = fv + rop(l022 +v4)+rop(l002 +v4)
      f4 = fv + rop(l021 +v4)+rop(l001 +v4)

      gradW_ny=(f2-f1)*tiy
      gradW_nz=(f4-f3)*tjz

      ! --- Gradients en T_j
      f1 =rop(il +v5)
      f2 =rop(ir +v5)

      gradT_ny=(f2-f1)*tiy

      !--- assemblage flux
      div = ( gradU_nx+gradV_ny+ gradW_nz)*cvisq
  
      f2 =  gradU_ny + gradV_nx
      f4 =  gradV_ny - div        
      f5 =  gradV_nz + gradW_ny 

#include  "FastS/Compute/mut_interface.for"

      fv     = -(    f2*tcy)*xmutvol
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     = -( 2.*f4*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     = -(    f5*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v4)+rop(il +v4))
      flu4   = flu4 + fv

      fv     = fv5*0.5  - (  gradT_ny*tcy )*xktvol
      flu5   = flu5 + fv

