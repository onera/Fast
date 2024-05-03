       ventx =  ventk(lven+v1ven)
       venty =  ventk(lven+v2ven)
       ventz =  ventk(lven+v3ven)*ck_vent

c       lven1 = lven +incj_ven              
c       lven2 = lven +inci_ven             
c       lven3 = lven2+incj_ven            

c       ventx = 0.25*( venti(lven  +v1ven)  !(i  , j,  k)
c     &               +venti(lven1 +v1ven)  !(i  , j+1,k)
c     &               +venti(lven2 +v1ven)  !(i+1, j  ,k)
c     &               +venti(lven3 +v1ven) )!(i+1, j  ,k)
c       venty = 0.25*( venti(lven  +v2ven)  
c     &               +venti(lven1 +v2ven)  
c     &               +venti(lven2 +v2ven)  
c     &               +venti(lven3 +v2ven) )
c       ventz = 0.25*( venti(lven  +v3ven)  
c     &               +venti(lven1 +v3ven) 
c     &               +venti(lven2 +v3ven)
c     &               +venti(lven3 +v3ven) )*ck_vent

