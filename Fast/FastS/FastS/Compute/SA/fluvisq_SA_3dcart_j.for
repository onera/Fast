      !--- Gradients  nutilde
      f1 =rop(il +v6)
      f2 =rop(ir +v6)

      gradT_ny=(f2-f1)*tiy

#include  "FastS/Compute/SA/mut_interface.for" 
      flu6 =flu6 - xmutvol*( gradT_ny*tcy )
