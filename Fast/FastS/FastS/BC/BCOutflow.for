      !l'etat exterieur est extrapole de l interieur sauf pour pression
      roe = rop(l1,1)
      rue = rop(l1,2)*roe
      rve = rop(l1,3)*roe
      rwe = rop(l1,4)*roe

      roe_inv = 1./state(1)

      ete = state(5) -0.5*( state(2)*state(2)
     &                +state(3)*state(3)
     &                +state(4)*state(4) )*roe_inv
     &          +0.5*(rue*rop(l1,2)+rve*rop(l1,3)+rwe*rop(l1,4) )

#include  "FastS/BC/non_reflection.for"

      rop(l,1) =  rop(l1,1)

      roinv   = 1./rop(l,1)

      if(qn.lt.0) then 
        rop(l,2)= rop(l1,2) - qn/(sn*sn)*tcx
        rop(l,3)= rop(l1,3) - qn/(sn*sn)*tcy
        rop(l,4)= rop(l1,4) - qn/(sn*sn)*tcz
      else 
        rop(l,2) =  rop(l1,2)
        rop(l,3) =  rop(l1,3) 
        rop(l,4) =  rop(l1,4) 
      endif

      rop(l,5) =  (qvar5*roinv - .5*( rop(l,2)*rop(l,2)
     &                               +rop(l,3)*rop(l,3)
     &                               +rop(l,4)*rop(l,4)) )*cvinv

