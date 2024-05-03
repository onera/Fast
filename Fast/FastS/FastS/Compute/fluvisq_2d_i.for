        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        lt200 = lt -  inci_mtr                    
        lt100 = lt +  inci_mtr                    
        lt010 = lt              +  incj_mtr
        lt210 = lt -  inci_mtr  +  incj_mtr 

        lvor  = lt
        lvol  = lt200   

        il   = l  -  inci
        ir   = l
        l220 = l  -  inci  -  incj
        l020 = l           -  incj
        l210 = l  -  inci  +  incj
        l010 = l           +  incj

CCCCCC
CCCCCC
CCCCCC   Facette I
CCCCCC
CCCCCC
      !Calcul des vecteurs surfaces I
      tix  = .5 * ( ti(lt + v1mtr) + ti(lt200+ v1mtr) )
      tiy  = .5 * ( ti(lt + v2mtr) + ti(lt200+ v2mtr) )
      tix1 = .5 * ( ti(lt + v1mtr) + ti(lt100+ v1mtr) )
      tiy1 = .5 * ( ti(lt + v2mtr) + ti(lt100+ v2mtr) )

      tjx  = .125 * ( tj(lt + v1mtr  ) + tj(lt200+ v1mtr) )
      tjy  = .125 * ( tj(lt + v2mtr  ) + tj(lt200+ v2mtr) )
      tjx1 = .125 * ( tj(lt010+ v1mtr) + tj(lt210+ v1mtr) )
      tjy1 = .125 * ( tj(lt010+ v2mtr) + tj(lt210+ v2mtr) )

      ! --- Gradients en U_i
      f1 = rop(il +v2)
      f2 = rop(ir +v2)
      f3 = rop(ir +v2)+rop(il +v2)+rop(l220 +v2)+rop(l020 +v2)
      f4 = rop(ir +v2)+rop(il +v2)+rop(l210 +v2)+rop(l010 +v2)

      gradU_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradU_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)

      !--- Gradients en v_i
      f1 = rop(il +v3)
      f2 = rop(ir +v3)
      f3 = rop(ir +v3)+rop(il +v3)+rop(l220 +v3)+rop(l020 +v3)
      f4 = rop(ir +v3)+rop(il +v3)+rop(l210 +v3)+rop(l010 +v3)

      gradV_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradV_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)

      !--- Gradients en T_i
      f1 = rop(il +v5)
      f2 = rop(ir +v5)
      f3 = rop(ir +v5)+rop(il +v5)+rop(l220 +v5)+rop(l020 +v5)
      f4 = rop(ir +v5)+rop(il +v5)+rop(l210 +v5)+rop(l010 +v5)

      gradT_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradT_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)

      !--- assemblage flux
      div = ( gradU_nx+gradV_ny)*cvisq
  
      f1 =  gradU_nx - div
      f2 =  gradU_ny + gradV_nx
      f4 =  gradV_ny - div        

      volinv = 1./(vol(lvor)+vol(lvol))
#include  "FastS/Compute/mut_interface.for"

      fv     = -(2.*f1*tcx +    f2*tcy)*xmutvol
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     = -(   f2*tcx + 2.*f4*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     = fv5*0.5  - ( gradT_nx*tcx + gradT_ny*tcy )*xktvol
      flu5   = flu5 + fv

