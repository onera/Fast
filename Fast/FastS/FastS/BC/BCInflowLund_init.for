      ijkplanrec= int(lund_param(2))
      clund     = lund_param(3)
      if( int(lund_param(1)).eq.0) clund =0.

      jamor     = int(lund_param(4))

      idirlund  = int(lund_param(5))  ! direction normal paroi

      if(jamor.eq.-1) jamor   = max(param_int(IJKV),param_int(IJKV+1),
     &                              param_int(IJKV+2))

      ci_amor = 0
      cj_amor = 0
      ck_amor = 0
      if (idirlund.eq.1) then
          ijk_amor = param_int(IJKV)
          ci_amor  = 1
      elseif(idirlund.eq.2) then
          ijk_amor = param_int(IJKV+1)
          cj_amor  = 1
      else
          ijk_amor = param_int(IJKV+2)
          ck_amor  = 1
      endif
