         flag_correct_flu=0

         !Condition limite paroi adiabatique 
         If(    bc_type.eq.BCWALL        .or.bc_type.eq.BCWALLVISCOUS
     &      .or.bc_type.eq.BCWALLVISCOUS_ISOT_FICH
     &      .or.bc_type.eq.BCWALLVISCOUS_TRANSITION) Then

            mobile_coef = 1.
            iptdata = param_int( param_int(PT_BC) + 1 + ndf + nb_bc)
            if (nbdata.ne.0) mobile_coef = param_real( iptdata )

            if(idir.le.2) THEN
             incijk     =1
             call bflwall(ndom, idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, ti,venti)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwall(ndom, idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tj,ventj)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwall(ndom, idir, mobile_coef, param_int(NEQ_K),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tk,ventk)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1

         !BCwallslip ou symmetry
         Elseif(bc_type.eq.BCWALLINVISCID.or.
     &          bc_type.eq.BCSYMMETRYPLANE) then 

            mobile_coef = 1.
            iptdata = param_int( param_int(PT_BC) + 1 + ndf + nb_bc)
            if (nbdata.ne.0) mobile_coef = param_real( iptdata )

            if(idir.le.2) THEN
             incijk     =1
             call bflwallslip(ndom,idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, ti,venti)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflwallslip(ndom,idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tj,ventj)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflwallslip(ndom,idir,mobile_coef, param_int(NEQ_K),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tk,ventk)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1


         !BCInj1 ou BCInflow    
         Elseif(bc_type.eq.BCINFLOW.or.bc_type.eq.BCINJ1) Then

            mobile_coef = 1.

            if(idir.le.2) THEN
             incijk     =1
             call bflinj1(ndom, idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, ti,venti)
            elseif(idir.le.4) THEN
             incijk     = param_int(NIJK)
             call bflinj1(ndom, idir, mobile_coef, param_int(NEQ_IJ),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tj,ventj)
            else
             incijk     = param_int(NIJK+1)*param_int(NIJK)
             call bflinj1(ndom, idir, mobile_coef, param_int(NEQ_K),
     &                    param_int, param_real, incijk, ind_CL,
     &                    rop,drodm, tk,ventk)
            endif

            !Mise a jour flag pour calcul correct des efforts
            flag_correct_flu=1

         Else
            continue
         Endif
