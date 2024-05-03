      !--- Gradients  nutilde
      f1 =rop(il +v6)
      f2 =rop(ir +v6)
      f5 =rop(ir +v6)+rop(il +v6)+rop(l220 +v6)+rop(l200 +v6)
      f6 =rop(ir +v6)+rop(il +v6)+rop(l120 +v6)+rop(l100 +v6)

      gradT_nx=(f2*tix1-f1*tix)+(f6*tkx1-f5*tkx)
      gradT_ny=(f2*tiy1-f1*tiy)+(f6*tky1-f5*tky)

#include  "FastS/Compute/SA/mut_interface.for"  
      flu6 =flu6 - xmutvol*( gradT_nx*tcx + gradT_ny*tcy )
