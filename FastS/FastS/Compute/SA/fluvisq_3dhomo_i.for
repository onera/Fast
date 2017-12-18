        !1=?+1; 
        !3=?+2 
        !2=?-1
        !4=?-2
        !l341  = (i+2,j-2,k+1)
        lt200 = lt -  inci_mtr                    
        lt100 = lt +  inci_mtr                    
        lt010 = lt              +  incj_mtr
        lt210 = lt -  inci_mtr  +  incj_mtr 

c        lvo   = lt
c        lvo200= lt200

        lvor  = lt
        lvol  = lt200

        il   = l  -  inci
        ir   = l

c        il = l  -  inci
        l220 = l  -  inci  -  incj
        l020 = l           -  incj
        l210 = l  -  inci  +  incj
        l010 = l           +  incj
        l202 = l  -  inci            -  inck
        l002 = l                     -  inck
        l201 = l  -  inci            +  inck
        l001 = l                     +  inck

CCCCCC
CCCCCC
CCCCCC   Facette I
CCCCCC
CCCCCC
      !Calcul des vecteurs surfaces I
      tix  = .5 * ( ti(lt    +v1mtr) + ti(lt200 +v1mtr) )
      tiy  = .5 * ( ti(lt    +v2mtr) + ti(lt200 +v2mtr) )
      tix1 = .5 * ( ti(lt    +v1mtr) + ti(lt100 +v1mtr) )
      tiy1 = .5 * ( ti(lt    +v2mtr) + ti(lt100 +v2mtr) )

      tjx  = .125 * ( tj(lt    +v1mtr) + tj(lt200 +v1mtr) )
      tjy  = .125 * ( tj(lt    +v2mtr) + tj(lt200 +v2mtr) )
      tjx1 = .125 * ( tj(lt010 +v1mtr) + tj(lt210 +v1mtr) )
      tjy1 = .125 * ( tj(lt010 +v2mtr) + tj(lt210 +v2mtr) )

      tkz  = .125 * ( tk(lt    +v1mtr) + tk(lt200 +v1mtr) )

      ! --- Gradients en U_i
      f1 = rop(il +v2)
      f2 = rop(ir +v2)
      f3 = rop(ir +v2)+rop(il +v2)+rop(l220 +v2)+rop(l020 +v2)
      f4 = rop(ir +v2)+rop(il +v2)+rop(l210 +v2)+rop(l010 +v2)
      f5 = rop(ir +v2)+rop(il +v2)+rop(l202 +v2)+rop(l002 +v2)
      f6 = rop(ir +v2)+rop(il +v2)+rop(l201 +v2)+rop(l001 +v2)

      gradU_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradU_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)
      gradU_nz=(f6-f5)*tkz

      !--- Gradients en v_i
      f1 = rop(il +v3)
      f2 = rop(ir +v3)
      f3 = rop(ir +v3)+rop(il +v3)+rop(l220 +v3)+rop(l020 +v3)
      f4 = rop(ir +v3)+rop(il +v3)+rop(l210 +v3)+rop(l010 +v3)
      f5 = rop(ir +v3)+rop(il +v3)+rop(l202 +v3)+rop(l002 +v3)
      f6 = rop(ir +v3)+rop(il +v3)+rop(l201 +v3)+rop(l001 +v3)

      gradV_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradV_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)
      gradV_nz=(f6-f5)*tkz

      !--- Gradients en W_i
      f1 = rop(il +v4)
      f2 = rop(ir +v4)
      f3 = rop(ir +v4)+rop(il +v4)+rop(l220 +v4)+rop(l020 +v4)
      f4 = rop(ir +v4)+rop(il +v4)+rop(l210 +v4)+rop(l010 +v4)
      f5 = rop(ir +v4)+rop(il +v4)+rop(l202 +v4)+rop(l002 +v4)
      f6 = rop(ir +v4)+rop(il +v4)+rop(l201 +v4)+rop(l001 +v4)

      gradW_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradW_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)
      gradW_nz=(f6-f5)*tkz

      !--- Gradients en T_i
      f1 = rop(il +v5)
      f2 = rop(ir +v5)
      f3 = rop(ir +v5)+rop(il +v5)+rop(l220 +v5)+rop(l020 +v5)
      f4 = rop(ir +v5)+rop(il +v5)+rop(l210 +v5)+rop(l010 +v5)

      gradT_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradT_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)

      !--- assemblage flux
      div = ( gradU_nx+gradV_ny+ gradW_nz)*cvisq
  
      f1 =  gradU_nx - div
      f2 =  gradU_ny + gradV_nx
      f3 =  gradU_nz + gradW_nx
      f4 =  gradV_ny - div        
      f5 =  gradV_nz + gradW_ny 

      volinv = 1./(vol(lvor)+vol(lvol))
#include  "FastS/Compute/mut_prandtltb_interface.for"

      fv     = -(2.*f1*tcx +    f2*tcy)*xmutvol
      fv5    =        fv * (rop(ir+v2)+rop(il +v2))
      flu2   = flu2 + fv

      fv     = -(   f2*tcx + 2.*f4*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v3)+rop(il +v3))
      flu3   = flu3 + fv

      fv     = -(   f3*tcx +    f5*tcy)*xmutvol
      fv5    = fv5  + fv * (rop(ir+v4)+rop(il +v4))
      flu4   = flu4 + fv

      fv     = fv5*0.5  - ( gradT_nx*tcx + gradT_ny*tcy )*xktvol
      flu5   = flu5 + fv
