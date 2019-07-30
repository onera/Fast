      ! on amorti les fluctu en dehors de la couche limite: si ijk > jamor defini dans muse
      amor  = 1
      iamor = ci_amor*i +cj_amor*j + ck_amor*k

      if (iamor.gt.jamor) then
        amor=abs(ijk_amor-iamor)/float(ijk_amor+1-jamor)
        amor=amor*amor
      endif

      !rop(l,1)=ro_inflow(li,1)+ clund*(rop(lr,1)-ro_planlund(li,1))*amor
      !rop(l,2)=ro_inflow(li,2)+ clund*(rop(lr,2)-ro_planlund(li,2))*amor
      rop(l,1)=ro_inflow(li,1)
      rop(l,2)=ro_inflow(li,2)
      rop(l,3)=ro_inflow(li,3)+ clund*(rop(lr,3)-ro_planlund(li,3))*amor
      rop(l,4)=ro_inflow(li,4)+ clund*(rop(lr,4)-ro_planlund(li,4))*amor

      rop(l,5) =  rop(l1,5)*rop(l1,1)/rop(l,1)
