        tcx  = tijk(lmtr,ic)*ci_mtr*snorm
        tcy  = tijk(lmtr,jc)*cj_mtr*snorm
        tcz  = tijk(lmtr,kc)*ck_mtr*snorm

        sn  = sqrt(tcx*tcx + tcy*tcy +tcz*tcz)

        roi = rop( l1 ,1)
        rui = rop( l1 ,2)*roi
        rvi = rop( l1 ,3)*roi
        rwi = rop( l1 ,4)*roi
        eti = rop( l1 ,5)*param_real(CVINF)*roi 
     &       +.5*(rui*rui+rvi*rvi+rwi*rwi)/roi

C Vitesse entrainement interieure

        qen =(  ventijk( ldp ,   1   )*tcx
     &        + ventijk( ldp ,   2   )*tcy
     &        + ventijk( ldp ,kc_vent)*tcz*ck_vent
     &       )*c_ale


c-----Define average of inside and outside conditions
c     as linearization state

        qref1 =.5*(roe+roi)
        qref2 =.5*(rue+rui)
        qref3 =.5*(rve+rvi)
        qref4 =.5*(rwe+rwi)
        qref5 =.5*(ete+eti)
c     
        
        r     = qref1
        ri    = 1./r
        u     = qref2*ri
        v     = qref3*ri
        w     = qref4*ri
        p     = gamm1*(qref5-.5*r*(u*u+v*v+w*w))
        c     = sqrt(param_real(GAMMA)*p*ri)
        ci    = 1./c
        s_1   = 1./sn
        qn    = tcx*u + tcy*v + tcz*w - qen

        if(sn.eq.0.)then
          tnx=1.
          tny=0.
        else
          tnx= tcx/sn
          tny= tcy/sn
        endif
c
c-----Compute eigenvalues from the average state
c
        eigv1=qn
        eigv2=qn
        eigv3=qn
        eigv4=qn+sn*c
        eigv5=qn-sn*c
c
c-----Transform outside conditions to characteristic variables
c
        qvar1=roe
        qvar2=rue
        qvar3=rve
        qvar4=rwe
        qvar5=ete
c
c-----Transforms from ordinary to characteristic variables
c
        svar1= qvar1
        svar2=(qvar2-u*qvar1)*ri
        svar3=(qvar3-v*qvar1)*ri
        svar4=(qvar4-w*qvar1)*ri
        svar5=gamm1*( .5*(u*u+v*v+w*w)*qvar1
     &               - u*qvar2-v*qvar3-w*qvar4
     &               + qvar5 )

        rvar1= svar1-ci**2*svar5
        rvar2= tcz*s_1*(tnx*svar2+tny*svar3) - svar4
        rvar3=         -tny*svar2+tnx*svar3
        rvar4= .5*s_1*( tcx*svar2+ tcy*svar3+ tcz*svar4 )+.5*ri*ci*svar5
        rvar5=-.5*s_1*( tcx*svar2+ tcy*svar3+ tcz*svar4 )+.5*ri*ci*svar5

        roext=rvar1
        ruext=rvar2
        rvext=rvar3
        rwext=rvar4
        etext=rvar5
c
c-----Transform inside conditions to characteristic variables
c
        qvar1=roi
        qvar2=rui
        qvar3=rvi
        qvar4=rwi
        qvar5=eti
c
c-----Transforms from ordinary to characteristic variables
c
        svar1= qvar1
        svar2=(qvar2-u*qvar1)*ri
        svar3=(qvar3-v*qvar1)*ri
        svar4=(qvar4-w*qvar1)*ri
        svar5=gamm1*( .5*(u*u+v*v+w*w)*qvar1
     .       - u*qvar2-v*qvar3-w*qvar4+qvar5 )
c
        rvar1= svar1-ci**2*svar5
        rvar2= tcz*s_1*(tnx*svar2+tny*svar3) - svar4
        rvar3=         -tny*svar2+tnx*svar3
        rvar4= .5*s_1*(tcx*svar2+ tcy*svar3+ tcz*svar4) +.5*ri*ci*svar5
        rvar5=-.5*s_1*(tcx*svar2+ tcy*svar3+ tcz*svar4) +.5*ri*ci*svar5
c

        roint=rvar1
        ruint=rvar2
        rvint=rvar3
        rwint=rvar4
        etint=rvar5
c
c-----Choose char. var. from inside or outside depending
c     on sign of corresponding eigenvalue
c
        rvar1=roint
        if(eigv1.lt.0.) rvar1=roext
        rvar2=ruint
        if(eigv2.lt.0.) rvar2=ruext
        rvar3=rvint
        if(eigv3.lt.0.) rvar3=rvext
        rvar4=rwint
        if(eigv4.lt.0.) rvar4=rwext
        rvar5=etint
        if(eigv5.lt.0.) rvar5=etext

c
c-----Transforms from characteristic to ordinary variables
c
        svar1= rvar1+r*ci*(rvar4+rvar5)
        svar2= tcz*s_1*tnx*rvar2-tny*rvar3 + tcx*s_1*(rvar4-rvar5)
        svar3= tcz*s_1*tny*rvar2+tnx*rvar3 + tcy*s_1*(rvar4-rvar5)
        svar4= -rvar2                      + tcz*s_1*(rvar4-rvar5)
        svar5= r*c*(rvar4+rvar5)
c
        qvar5=.5*(u*u+v*v+w*w)*svar1    + r*(u*svar2+v*svar3+w*svar4)
     &                                  +gamm1_1*svar5

        !roinv = 1./svar1
        !u  = u + r*svar2*roinv
        !v  = v + r*svar3*roinv
        !w  = w + r*svar4*roinv

       !rop(l,1) =  svar1
       !rop(l,2) =  u 
       ! rop(l,3) =  v 
       ! rop(l,4) =  w 
       ! rop(l,5) =  (qvar5*roinv - .5*(u*u+v*v+w*w))*cvinv

