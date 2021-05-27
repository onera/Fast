         flag_correct_flu=0
          
         !Wall model boundary condition
         If(bc_type.eq.BCWALLMODEL.or.bc_type.eq.BCWALLEXCHANGE) then
            
            nbdata = param_int(pt_bc + BC_NBDATA)

            if(idir.le.2) THEN
             incijk     =1
             call bflwallmodel(ndom,ithread, idir, nbdata, nitcfg,
     &                    param_int(NEQ_IJ), param_real(iptdata),
     &                    param_int, param_real,
     &                    incijk, ind_CL,rop,drodm, ti,venti,xmut,x,y,z)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwallmodel(ndom,ithread, idir, nbdata, nitcfg,
     &                    param_int(NEQ_IJ), param_real(iptdata),
     &                    param_int, param_real,
     &                    incijk, ind_CL, rop,drodm,tj,ventj,xmut,x,y,z)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwallmodel(ndom,ithread, idir, nbdata, nitcfg,
     &                    param_int(NEQ_K), param_real(iptdata),
     &                    param_int, param_real,
     &                    incijk, ind_CL, rop,drodm,tk,ventk,xmut,x,y,z)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1
         Else
            continue
         Endif
