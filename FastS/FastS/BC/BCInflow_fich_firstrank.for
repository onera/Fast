      !l'etat  impose, sauf pression extrapolé de l'interieur
       !roi0    = pa(li)*param_real(GAMMA)*gamm1_1/ha(li)
       !roe     = roi0*(1+ gamm1*0.5*mach*mach)**(-gamm1_1)

       !modif ivan
       roe     = pa(li)
       roe_inv = 1./roe       

       ru  = sqrt( 2.*(pa(li)-pref)*roe)

       !modif ivan
       !rue=ru*d0x(li)
       !rve=ru*d0y(li)
       !rwe=ru*d0z(li)
       !!
       !! dox=velocityX,...., Pa=Density

       !!
       rue=d0x(li)*pa(li)
       rve=d0y(li)*pa(li)
       rwe=d0z(li)*pa(li)
       !rve=0.
       !rwe=0.

 
      ete = rop(l1,5)* rop(l1,1)*param_real(CVINF)
     &     +0.5*( rue*rue + rve*rve + rwe*rwe) *roe_inv

#include  "FastS/BC/non_reflection.for"

      roinv = 1./svar1

      u  = u + r*svar2*roinv
      v  = v + r*svar3*roinv
      w  = w + r*svar4*roinv


      !rop(l,1) =  svar1
      !rop(l,2) =  u 
      !rop(l,3) =  v 
      !rop(l,4) =  w 
      rop(l,1) =  roe
      rop(l,2) =  d0x(li) 
      rop(l,3) =  d0y(li)
      rop(l,4) =  d0z(li) 
      !rop(l,5) =  (qvar5*roinv - .5*(u*u+v*v+w*w))*cvinv
      rop(l,5) =  rop(l1,5)*rop(l1,1)*roe_inv 
