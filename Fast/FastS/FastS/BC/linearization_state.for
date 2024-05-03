C ... Implementation elsA       
        tcx  = tijk(lmtr,ic)*ci_mtr*snorm
        tcy  = tijk(lmtr,jc)*cj_mtr*snorm
        tcz  = tijk(lmtr,kc)*ck_mtr*snorm

        sn  = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)

        s_1   = 1./sn
        
        if(sn.eq.0.)then
          tnx=1.
          tny=0.
          tnz=0.
        else
          tnx= tcx*s_1
          tny= tcy*s_1
          tnz= tcz*s_1
        endif

C ... Inner state 
        roi = rop( l1 ,1)
        ui  = rop( l1 ,2)
        vi  = rop( l1 ,3)
        wi  = rop( l1 ,4)
        Ti  = rop( l1 ,5)
        pi  = roi*rgp*Ti
        
C ... linearization state = inner state (elsA-like)
        ro0 = rop( l1 ,1)
        u0  = rop( l1 ,2)
        v0  = rop( l1 ,3)
        w0  = rop( l1 ,4)
        T0  = rop( l1 ,5)
        p0  = ro0*rgp*T0
