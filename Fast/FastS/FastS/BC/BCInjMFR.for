CC    IN : rop
CC    OUT: roi,ui,vi,wi,Ti,pi (etat intern)
CC    OUT: ro0,u0,v0,w0,T0,pi (etat extern)
CC    OUT: tcx,tcy,tcz, tnx,tny,tnz
#include  "FastS/BC/linearization_state.for"
CC    In : tnx,tny,tnz, c_ale, ventijk
CC    OUT: qen (vitesse entrainement)
#include  "FastS/BC/qn_ale.for"
C ... //////////////////////////////////////////////////////////////////////////
C ... Polynomial solving
C qp(size_data)    Mass flow per area unit (B.C) defined on the interfaces
C ha(size_data)    Total enthalpy
C d0x(size_data)   x-component of the flow direction (B.C)
C d0y(size_data)   y-component of the flow direction (B.C)
C d0z(size_data)   z-component of the flow direction (B.C)
C nuext(size_data) turbulent viscosity
C ... //////////////////////////////////////////////////////////////////////////
C ...  Implementation elsA
        roc0 = SQRT(ro0*gam*p0)
        usd0n = 1/(d0x(li)*tnx + d0y(li)*tny + d0z(li)*tnz)
        usd0n2 = usd0n*usd0n
        wn = u0*tnx + v0*tny + w0*tnz - qen
        IF (ABS(wn) < 1.e-6) THEN
c ... Pressure newton computation -> seulement pour "amorcer"
          residug  = 1.e+20
          nitnwt   = 0
          mach = SQRT((u0*u0)+(v0*v0)+(w0*w0))
          mach = mach*ro0/roc0
          p  = p0*(1+0.2*(mach*mach))**gam2
          ri = pi + roc0*wn
          DO WHILE (residug .GT. tolnewton .AND. nitnwt .LT. newtonmax)
            nitnwt  = nitnwt+1
            residug = 0.
            b = max(0.2,1. - ((wn+qen)**2)*usd0n2/(2.*ha(li)))
            pg = p*b**gam2
            rog = gam2*pg/(ha(li)*b)
            f = pg + roc0*wn - ri
            df = roc0 - rog*(wn+qen)*usd0n2
            dwng= -f/df
            wn = wn + dwng
            residug = max(residug,ABS(dwng/wn))
          ENDDO
          qn = wn+qen
          b = max(0.2,1. - (qn**2)*usd0n2/(2.*ha(li)))
          pg = p*b**gam2
          rog = gam2*pg/(ha(li)*b)
          Tg = pg/(rog*rgp)
          w = qn*usd0n
          rop(l,1) = rog          ! ro
          rop(l,2) = w*d0x(li)    ! u
          rop(l,3) = w*d0y(li)    ! v
          rop(l,4) = w*d0z(li)    ! w
          rop(l,5) = Tg           ! T
        ELSE
c   ... Solving the unknown "qn" (1-state normal velocity)
c   ... Polynomial expression : ax**2 + bx + c = 0
          coefa = 0.5*qp(li)*usd0n2 / gam2 + roc0              ! 1/2*q*1/(d.n)**2*(gam-1)/gam + ro*c
          coefb = p0 + roc0*(u0*tnx + v0*tny + w0*tnz - qen)  ! p + ro*c*Un (normale sortante)
          coefc = -ha(li)*qp(li)/gam2                         ! - ha*q*(gam-1)/gam
          delta = coefb*coefb - 4.*coefa*coefc
          qn = ( -coefb + SQRT(delta))/(2.*coefa) + qen
          w = qn*ABS(usd0n)
          rop(l,1) = qp(li)/qn                                   ! ro
          rop(l,2) = w*d0x(li)                             ! u
          rop(l,3) = w*d0y(li)                             ! v
          rop(l,4) = w*d0z(li)                             ! w
          rop(l,5) = (ha(li)-0.5*w*w)/(param_real(CVINF)*gam)      ! T -> (ha-1/2*U**2)/(Cp)
        ENDIF


C ... //////////////////////////////////////////////////////////////////////////
C ... //////////////////////////////////////////////////////////////////////////
