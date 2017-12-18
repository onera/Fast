        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        lt100 = lt +  inci_mtr                    
        lt010 = lt            +  incj_mtr                      
        lt002 = lt                           !pas constant en z

        lvor  = lt
        lvol  = lt002

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
      tiz  =  tk(lt    +v1mtr) 

      tjx  =  0.25*ti(lt    +v1mtr)
      tjy  =  0.25*ti(lt    +v2mtr)
      tjx1 =  0.25*ti(lt100 +v1mtr)
      tjy1 =  0.25*ti(lt100 +v2mtr)

      tkx  =  0.25*tj(lt    +v1mtr)
      tky  =  0.25*tj(lt    +v2mtr)
      tkx1 =  0.25*tj(lt010 +v1mtr)
      tky1 =  0.25*tj(lt010 +v2mtr)

      ! --- Gradients en U_k
      f1 = rop(il +v2)
      f2 = rop(ir +v2)
      f3 = rop(ir +v2)+rop(il +v2)+rop(l202 +v2)+rop(l200 +v2)
      f4 = rop(ir +v2)+rop(il +v2)+rop(l102 +v2)+rop(l100 +v2)
      f5 = rop(ir +v2)+rop(il +v2)+rop(l022 +v2)+rop(l020 +v2)
      f6 = rop(ir +v2)+rop(il +v2)+rop(l012 +v2)+rop(l010 +v2)

      gradU_nx=((f4*tjx1-f3*tjx)+(f6*tkx1-f5*tkx))
      gradU_nz=(f2-f1)*tiz

      ! --- Gradients en V_k
      f1 = rop(il +v3)
      f2 = rop(ir +v3)
      f3 = rop(ir +v3)+rop(il +v3)+rop(l202 +v3)+rop(l200 +v3)
      f4 = rop(ir +v3)+rop(il +v3)+rop(l102 +v3)+rop(l100 +v3)
      f5 = rop(ir +v3)+rop(il +v3)+rop(l022 +v3)+rop(l020 +v3)
      f6 = rop(ir +v3)+rop(il +v3)+rop(l012 +v3)+rop(l010 +v3)

      gradV_ny=((f4*tjy1-f3*tjy)+(f6*tky1-f5*tky))
      gradV_nz=(f2-f1)*tiz

      ! --- Gradients en W_k
      f1 = rop(il +v4)
      f2 = rop(ir +v4)
      f3 = rop(ir +v4)+rop(il +v4)+rop(l202 +v4)+rop(l200 +v4)
      f4 = rop(ir +v4)+rop(il +v4)+rop(l102 +v4)+rop(l100 +v4)
      f5 = rop(ir +v4)+rop(il +v4)+rop(l022 +v4)+rop(l020 +v4)
      f6 = rop(ir +v4)+rop(il +v4)+rop(l012 +v4)+rop(l010 +v4)

      gradW_nx=((f4*tjx1-f3*tjx)+(f6*tkx1-f5*tkx))
      gradW_ny=((f4*tjy1-f3*tjy)+(f6*tky1-f5*tky))
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

      volinv = 1./(vol(lvor)+vol(lvol))
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
