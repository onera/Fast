      li1  = indcg(ind_loop(1)  , 1, 1)
      li2  = indcg(ind_loop(2)+1, 1, 1)

      xlong=sqrt( (x(li1)-x(li2))**2+(y(li1)-y(li2))**2)

      ! Funk
      !zlong_transi   =  0.
      !freq_transi    =  80000.00
      !lambdaz_transi =  1.
      !ampli_transi   =  0.0005

      !Guillaume sd7003
      zlong_transi   =  0.
      freq_transi    =  2.3
      lambdaz_transi =  2.
      ampli_transi   =  0.0005

      if(zlong_transi.eq.0.0) then
        zlong=abs( z( indcg(ind_loop(1),1,       1              ) )
     &            -z( indcg(ind_loop(1),1, param_int(IJKV+2)+1) ) )
      else
        zlong=zlong_transi
      endif

      x0   = x(indcg((ind_loop(1)+ind_loop(2))/2,1,1))
      y0   = y(indcg((ind_loop(1)+ind_loop(2))/2,1,1))

      omega  = 2*pi*freq_transi
      lambda =2*pi*lambdaz_transi/zlong
      !write(*,*)'long transi',xlong,x0,'freq et long onde',omega,lambda

      xlong=1./xlong/xlong


c          !! mise a jour du random a la premiere ssiter uniquememnt
c          if(nstep.eq.1) then
c            do k = ind_loop(5), ind_loop(6)
c            do i = ind_loop(1), ind_loop(2)
c                li  = indbci(i,  k )
c                call random_number(random(li))
c            enddo
c            enddo
c          endif

#if __DEBUG__
            ! encadre le depassement intentionnel
            if (c_ale > 0) then
#endif
            ventx      = ventijk(ldp ,1)
            venty      = ventijk(ldp ,2)
            ventz      = ventijk(ldp ,kc_vent)*ck_vent
#if __DEBUG__
            else
            ventx = 0.
            venty = 0.
            ventz = 0.
            endif
#endif
            tcx  = tijk(lmtr,ic)*ci_mtr
            tcy  = tijk(lmtr,jc)*cj_mtr
            tcz  = tijk(lmtr,kc)*ck_mtr

          !on calcule l'inverse de la surface de la facette
          s_1  = c_pertu/sqrt(tcx*tcx + tcy*tcy +tcz*tcz)


          !! mise a jour du random a la premiere ssiter uniquememnt
          !rnd = random(li) 
          rnd = 0.
          
          !pertu random
          ! pertu sinuoidale temps et z + gaussien (x,y)
          rnd = (rnd-0.5)*ampli_transi*0.05             
     &       +2.*ampli_transi*sin(omega*param_real(TEMPS))
     &            *sin(lambda*z(ldx))      
     &            *exp(-10.*xlong*( (x(ldx)-x0)*(x(ldx)-x0)
     &                         +(y(ldx)-y0)*(y(ldx)-y0)))

            u  =  rop(ldjr,2) - tcx*rnd*s_1
            v  =  rop(ldjr,3) - tcy*rnd*s_1
            w  =  rop(ldjr,4) -     rnd*c_pertu

            rop(l,1) = rop(ldjr,1)
            rop(l,2) = 2.*ventx*c_ale - u
            rop(l,3) = 2.*venty*c_ale - v
            rop(l,4) = 2.*ventz*c_ale - w

            rop(l,5) = rop(ldjr,5)
