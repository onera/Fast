      qp1 = rop(nm +v1) 
      qp2 = rop(nm +v2) 
      qp3 = rop(nm +v3)
      qp4 = rop(nm +v4)
      qm1 = rop( l +v1) 
      qm2 = rop( l +v2) 
      qm3 = rop( l +v3)
      qm4 = rop( l +v4)
      
      !moyenne roe pour o1  
      qp1 = sqrt( qp1 )
      qm1 = sqrt( qm1 )
      r_1 = 1./(qp1 + qm1)
  
      u   = r_1*( qp1*qp2 + qm1*qm2)
      v   = r_1*( qp1*qp3 + qm1*qm3)
      w   = r_1*( qp1*qp4 + qm1*qm4)

      !ro_moy     = sqrt( rop(nm +v1)*rop( l +v1) )
 
      qn   = u*tcx+v*tcy+w*tcz -qen
      tdu  =  max(c1/s_1, abs(qn) )
      ! s_1 = 1./Surf
  
      qn2=(rop(l+v2)*tcx +rop(l+v3)*tcy +rop(l+v4)*tcz) -qen
      qm6=rop(l+v6)*rop(l+v1)*qn2
     
      qn1=(rop(nm+v2)*tcx +rop(nm+v3)*tcy +rop(nm+v4)*tcz) -qen
      qp6=rop(nm+v6)*rop(nm+v1)*qn1

      flu6 = 0.5*(qm6 + qp6 - tdu*( rop(l +v6)*rop(l +v1)
     &                             -rop(nm+v6)*rop(nm+v1)
     &                            )
     &          )

c             if(ndom.eq.1.or.ndom.eq.4) then
c               test1 =isnan(flu6)
c               if (test1) then
c                    write(*,*) 'flu6',l-lij+1,j,k
c                    stop
c               endif
c              endif

c      qm6  = qm*r2*qn2
c      qp6  = qp*r1*qn1

c      flu6 = 0.5*(qm6 + qp6 +sign(1.,flu1)*(qp6-qm6))
