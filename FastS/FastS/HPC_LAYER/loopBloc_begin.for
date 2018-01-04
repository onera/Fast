      lth  = 0
      test(3) = mod(ijkv_thread(3), skip(3))
      do kGbloc= 1, ijkv_thread(3), skip(3)

        shift(3)= skip(3)-1
        if(test(3).ne.0.and.kGbloc.eq.ijkv_thread(3))shift(3)=shift(3)-1

      test(2) = mod(ijkv_thread(2), skip(2))
      do jGbloc= 1, ijkv_thread(2), skip(2)

        shift(2)= skip(2)-1
        if(test(2).ne.0.and.jGbloc.eq.ijkv_thread(2))shift(2)=shift(2)-1

      test(1) = mod(ijkv_thread(1), skip(1))
      do iGbloc= 1, ijkv_thread(1), skip(1)

       !!
       !!
       !!si pas assez de travail pour tous les threads, on skip....
       !!
       !!
        if(ithread_sock.gt.thread_parsock_actif) goto 9999

        shift(1)= skip(1)-1
        if(test(1).ne.0.and.iGbloc.eq.ijkv_thread(1))shift(1)=shift(1)-1

#if CHECK_SPLIT > 1
         inc11=0
         inc22=0
#endif 
         do kp=0,shift(3)

          kbloc = kGbloc + kp
          sens(3)=1
          if(kp.ne.0)  sens(3)=-1

           do jp=0,shift(2)

            jbloc = jGbloc + jp
            sens(2)=1
            if(jp.ne.0)  sens(2)=-1

            do ip=0,shift(1)

            ibloc = iGbloc + ip
            sens(1)=1
            if(ip.ne.0)  sens(1)=-1

     
        !modifier ithred et thread_pos

cc
         thread_pos_tmp(1)= (1-2*ip)*thread_pos(1) + ip*(topo_s(1)+1) 
         thread_pos_tmp(2)= (1-2*jp)*thread_pos(2) + jp*(topo_s(2)+1) 
         thread_pos_tmp(3)= (1-2*kp)*thread_pos(3) + kp*(topo_s(3)+1) 

         ithread_sock = thread_pos_tmp(1)
     &                 +(thread_pos_tmp(2)-1)*topo_s(1)  
     &                 +(thread_pos_tmp(3)-1)*topo_s(1)*topo_s(2)  


        lth    = lth +1
        lwait  = 0
        lgo    = 0
        !if(ip.eq.shift(1).and.jp.eq.shift(2).and.kp.eq.shift(3)) lgo  =2
        !if(ip+jp+kp.eq.0.and.lth.ne.1)                           lwait=2
        if(ip.eq.shift(1).and.jp.eq.shift(2).and.kp.eq.shift(3)) lgo  =1
        if(ip+jp+kp.eq.0.and.lth.ne.1)                           lwait=1
        
        if(Nbre_socket.eq.1) then
          lgo   = -1
          lwait = -1
        endif
