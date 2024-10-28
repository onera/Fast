#ifndef E_SCALAR_COMPUTER
#include   "FastS/Compute/loopGpu_begin.for"
#else
           do k = ind_loop(5), ind_loop(6)
           do j = ind_loop(3), ind_loop(4)
           do i = ind_loop(1), ind_loop(2)

             l  = inddm(i,j,k)
             lt = indmtr(i,j,k)
#endif 
             lvo= lt
