      !--- Gradients  nutilde
      f1 = rop(il +v6)
      f2 = rop(ir +v6)
      f3 = rop(ir +v6)+rop(il +v6)+rop(l202 +v6)+rop(l200 +v6)
      f4 = rop(ir +v6)+rop(il +v6)+rop(l102 +v6)+rop(l100 +v6)
      f5 = rop(ir +v6)+rop(il +v6)+rop(l022 +v6)+rop(l020 +v6)
      f6 = rop(ir +v6)+rop(il +v6)+rop(l012 +v6)+rop(l010 +v6)

      gradT_nx=(f2*tix1-f1*tix)+ (f4*tjx1-f3*tjx)+(f6*tkx1-f5*tkx)
      gradT_ny=(f2*tiy1-f1*tiy)+ (f4*tjy1-f3*tjy)+(f6*tky1-f5*tky)
      gradT_nz=(f2*tiz1-f1*tiz)+ (f4*tjz1-f3*tjz)+(f6*tkz1-f5*tkz)

#include  "FastS/Compute/SA/assemble_tauij_3dfull.for"
