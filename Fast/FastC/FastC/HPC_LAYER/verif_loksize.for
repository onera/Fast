        !! on dimsension les verrous sur le param_int( IO_THREAD)mbre d'element par face
        !dans les directions de decoupage
        !on le fait au premier bloc seulement pour garantir que la forme
        !du tableau ne change. On prend une marge 2*2*2=8 pour garantir
        !une taille suffisante pour le dernier bloc

c        ipt_lok_sock = Nbre_thread_actif + 1
c        ipt_lok      = ipt_lok_sock  + lok_shap_sock(4)*lok_shap_sock(3)
c     &                                               *Nbre_socket
        ipt_lok =1

        lok_shap(1)  = 2* min( thread_topology(1)-1, 1 )
        lok_shap(2)  = 2* min( thread_topology(2)-1, 1 )
        lok_shap(3)  = 2* min( thread_topology(3)-1, 1 )
        
        size_max = ijkv_sdm(2)*ijkv_sdm(3)*lok_shap(1)
        size_loc = ijkv_sdm(1)*ijkv_sdm(3)*lok_shap(2)
        size_max = max(size_max,size_loc)
        size_loc = ijkv_sdm(1)*ijkv_sdm(2)*lok_shap(3)
        size_max = max(size_max,size_loc)

        !!socket=1, pas besoin de marge
        neq_lok    = 1
        lok_shap(3)= max(1, lok_shap(1)+lok_shap(2)+lok_shap(3))
        lok_shap(4)= mx_synchro/lok_shap(3)

        if(size_max .gt.lok_shap(4))  lerr=.true.
         
        if( lerr  ) then
          !$OMP   SINGLE
            write(*,*)'------'
            write(*,*)'error msg'
            write(*,*)'------'
            write(*,'(a,i5)')'resize MX_SYNCHRO. Present value=',
     & mx_synchro
            write(*,*)'Value must be at least larger than  :',
     & lok_shap(3)*size_max
            write(*,*)'Just after modules import of userscript.py,
     & add the following python command:'
            write(*,'(a)')'#'
            write(*,'(a)')'#'
            write(*,'(a,i5)')'import FastC.PyTree as FastC'
            write(*,'(a,i5)')'FastC.MX_SYNCHRO=',
     & lok_shap(3)*size_max+1
            write(*,*)'------'
            write(*,*)' End error msg'
            write(*,*)'------'
            call error('navier_stokes_struct$',70,1)
          !$OMP   END SINGLE
        endif

 
