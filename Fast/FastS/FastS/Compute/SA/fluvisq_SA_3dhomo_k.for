      !--- Gradients  nutilde
      f1 =rop(il +v6)
      f2 =rop(ir +v6)

      gradT_nz=(f2-f1)*tiz

#include  "FastS/Compute/SA/mut_interface.for"  
      flu6 =flu6 - xmutvol*gradT_nz*tcz
