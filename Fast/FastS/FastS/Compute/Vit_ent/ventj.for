       ventx =  ventj(lven+v1ven)
       venty =  ventj(lven+v2ven)
       ventz =  ventj(lven+v3ven)*ck_vent

c       lven1 = lven +inci_ven              
c       lven2 = lven +inck_ven             
c       lven3 = lven2+inci_ven            

c       ventx = 0.25*( venti(lven  +v1ven)  !(i  , j, k  )
c     &               +venti(lven1 +v1ven)  !(i+1, j, k  )
c     &               +venti(lven2 +v1ven)  !(i  , j, k+1)
c     &               +venti(lven3 +v1ven) )!(i  , j, k+1)
c       venty = 0.25*( venti(lven  +v2ven)  
c     &               +venti(lven1 +v2ven)  
c     &               +venti(lven2 +v2ven)  
c     &               +venti(lven3 +v2ven) )
c       ventz = 0.25*( venti(lven  +v3ven)  
c     &               +venti(lven1 +v3ven) 
c     &               +venti(lven2 +v3ven)
c     &               +venti(lven3 +v3ven) )*ck_vent

