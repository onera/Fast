              vis   = 1e30
              conv  = 1e30
              if (isSound .eq. 0) then
                 c=0
              endif
              if (isconv .eq. 1) then
                 xinvspc = 1./(sp+c)
                 conv    = 2.*vol(lvo) * surf * xinvspc
                 comp    = conv
              endif
              if (isvisc .eq. 1) then
                 xinvkt = gam1/xmut(l)
                 vis    = 2.*r*vol(lvo)*vol(lvo)*surf*surf*xinvkt
                 comp   = min(vis,conv)
              endif             
              cfloc = 1.0/comp
              cfl(1) = max(cfl(1),cfloc)
              cfl(2) = min(cfl(2),cfloc)
              cfl(3) = cfl(3) + cfloc
              ipt_cfl(l)=cfloc
