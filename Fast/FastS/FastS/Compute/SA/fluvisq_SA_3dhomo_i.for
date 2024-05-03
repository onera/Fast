      !--- Gradients  nutilde
      f1 = rop(il +v6)
      f2 = rop(ir +v6)
      f3 = rop(ir +v6)+rop(il +v6)+rop(l220 +v6)+rop(l020 +v6)
      f4 = rop(ir +v6)+rop(il +v6)+rop(l210 +v6)+rop(l010 +v6)

      gradT_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)
      gradT_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)

#include  "FastS/Compute/SA/mut_interface.for" 
      flu6 =flu6 - xmutvol*(  gradT_nx*tcx + gradT_ny*tcy  ) 
