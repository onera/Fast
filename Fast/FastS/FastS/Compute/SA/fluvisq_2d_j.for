        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        lt100 = lt +  inci_mtr                    
        lt010 = lt                +  incj_mtr   
        lt020 = lt                -  incj_mtr   
        lt120 = lt +  inci_mtr    -  incj_mtr                 

        lvor  = lt
        lvol  = lt020

        il    = l             -  incj
        ir    = l  
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
      tix  = .5 * ( tj(lt + v1mtr) + tj(lt020 + v1mtr) )
      tiy  = .5 * ( tj(lt + v2mtr) + tj(lt020 + v2mtr) )
      tix1 = .5 * ( tj(lt + v1mtr) + tj(lt010 + v1mtr) )
      tiy1 = .5 * ( tj(lt + v2mtr) + tj(lt010 + v2mtr) )

      tkx  = .125 * ( ti(lt    + v1mtr)+ ti(lt020 + v1mtr) ) !(i,j,k) et (i,j-1,k)
      tky  = .125 * ( ti(lt    + v2mtr)+ ti(lt020 + v2mtr) )
      tkx1 = .125 * ( ti(lt100 + v1mtr)+ ti(lt120 + v1mtr) ) !(i+1,j,k) et (i+1,j-1,k)
      tky1 = .125 * ( ti(lt100 + v2mtr)+ ti(lt120 + v2mtr) )

      ! --- Gradients en U_j (i,j,k), (i,j-1,k), (i,j-1,k-1), (i,j,k-1)
      f1 =rop(il +v2)
      f2 =rop(ir +v2)
      f5 =rop(ir +v2)+rop(il +v2)+rop(l220 +v2)+rop(l200 +v2)
      f6 =rop(ir +v2)+rop(il +v2)+rop(l120 +v2)+rop(l100 +v2)

      gradU_nx=(f2*tix1-f1*tix)+(f6*tkx1-f5*tkx)
      gradU_ny=(f2*tiy1-f1*tiy)+(f6*tky1-f5*tky)

      ! --- Gradients en V_j
      f1 =rop(il +v3)
      f2 =rop(ir +v3)
      f5 =rop(ir +v3)+rop(il +v3)+rop(l220 +v3)+rop(l200 +v3)
      f6 =rop(ir +v3)+rop(il +v3)+rop(l120 +v3)+rop(l100 +v3)

      gradV_nx=(f2*tix1-f1*tix)+(f6*tkx1-f5*tkx)
      gradV_ny=(f2*tiy1-f1*tiy)+(f6*tky1-f5*tky)

      ! --- Gradients en T_j
      f1 =rop(il +v5)
      f2 =rop(ir +v5)
      f5 =rop(ir +v5)+rop(il +v5)+rop(l220 +v5)+rop(l200 +v5)
      f6 =rop(ir +v5)+rop(il +v5)+rop(l120 +v5)+rop(l100 +v5)

      gradT_nx=(f2*tix1-f1*tix)+(f6*tkx1-f5*tkx)
      gradT_ny=(f2*tiy1-f1*tiy)+(f6*tky1-f5*tky)

      !--- assemblage flux
      div = ( gradU_nx+gradV_ny)*cvisq
  
      f1 =  gradU_nx - div
      f2 =  gradU_ny + gradV_nx
      f4 =  gradV_ny - div        

      volinv = 1./(vol(lvor)+vol(lvol))
#include  "FastS/Compute/mut_prandtltb_interface.for"

      fv     = -(2.*f1*tcx +    f2*tcy)*xmutvol
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     = -(   f2*tcx + 2.*f4*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     = fv5*0.5  - ( gradT_nx*tcx + gradT_ny*tcy )*xktvol
      flu5   = flu5 + fv
