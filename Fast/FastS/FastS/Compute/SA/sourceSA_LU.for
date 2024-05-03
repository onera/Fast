               !pour la Jacobienne (linearisation du terme visqueux)
               tsourceb  = prod-destruc
               amutild   = anutild*rop(l,1)
               amutild   = max(amutild,0.0000000000000001)      
               tsourcenu = vol(lvo)*tsourceb/amutild

               coe(l, icoe_pos) = coe(l,5) - min(tsourcenu,0.)*coe(l,1)
