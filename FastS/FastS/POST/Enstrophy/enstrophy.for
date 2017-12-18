           u3      = dvardc(lg+ Vdvdx)- dvardc(lg+ Vdudx+dim_grad)    !dvdx-dudy
           enst(l) = u3*u3

           u3      = dvardc(lg+ Vdudx+2*dim_grad)- dvardc(lg+ Vdwdx)  !dudz-dwdx
           enst(l) = enst(l)+ u3*u3

           u3      = dvardc(lg+ Vdwdx+  dim_grad)
     &              -dvardc(lg+ Vdvdx+2*dim_grad)  !dwdy-dvdz
           enst(l) = enst(l)+ u3*u3
           
           enst(l)= 0.5*rop(l+ v1)*enst(l)*volinv 

