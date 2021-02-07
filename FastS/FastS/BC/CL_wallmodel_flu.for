         flag_correct_flu=0
          
         !Wall model boundary condition
         If(bc_type.eq.BCWALLMODEL) then
            
            sampling = 1   ! Wall model input from first cell
            if(idir.le.2) THEN
             incijk     =1
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_IJ), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, ti,venti,xmut,
     &                    sampling)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_IJ), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, tj,ventj,xmut,
     &                    sampling)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_K), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, tk,ventk,xmut,
     &                    sampling)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1

         elseif(bc_type.eq.BCWALLEXCHANGE) then
            
            sampling = 4   ! Wall model input from fourth cell
            if(idir.le.2) THEN
             incijk     =1
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_IJ), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, ti,venti,xmut,
     &                    sampling)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_IJ), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, tj,ventj,xmut,
     &                    sampling)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwallmodel(ndom, idir, iptparam, nitcfg,
     &                    param_int(NEQ_K), param_int, param_real,
     &                    incijk, ind_CL, rop,drodm, tk,ventk,xmut,
     &                    sampling)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1
         Else
            continue
         Endif
