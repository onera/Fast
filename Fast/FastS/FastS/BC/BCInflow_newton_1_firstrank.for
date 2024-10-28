CC    IN : rop
CC    OUT: roi,ui,vi,wi,Ti,pi (etat intern)
CC    OUT: ro0,u0,v0,w0,T0,pi (etat extern)
CC    OUT: tcx,tcy,tcz, tnx,tny,tnz
#include  "FastS/BC/linearization_state.for"
CC    In : tnx,tny,tnz, c_ale, ventijk
CC    OUT: qen (vitesse entrainememnt)
#include  "FastS/BC/qn_ale.for"

C. ... Non-linear equation for the normal velocity => Newton
C ...  Implementation elsA 
        vni  = ui*tnx + vi*tny + wi*tnz
        
        roc0 = SQRT(ro0*gam*p0) ! rho*c
        
        !sans dimension: produit scalaire direction vitesse . normale
        usd0n  = 1./(d0x(li)*tnx + d0y(li)*tny + d0z(li)*tnz)
        usd0n2 = usd0n**2
        
C ...   Inner caracteristic variable
C ...   Relative normal velocity        
        wni       = vni - qen
        ri   = pi + roc0*wni
        
C ...   Newton Initialization for the relative normal velocity 
        wn0       = u0*tnx + v0*tny + w0*tnz  - qen       
        wng  = wn0
