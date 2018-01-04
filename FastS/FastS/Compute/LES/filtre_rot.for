             !plan k-1
             l0  = l  -inci -incj - inck  + var
             l1  = l0 +inci 
             l2  = l1 +inci

             l3  = l0       +incj 
             l4  = l3 +inci 
             l5  = l4 +inci 

             l6  = l3       +incj
             l7  = l4       +incj 
             l8  = l5       +incj 
             
             fv =  c111*( rot(l0) + rot(l2) + rot(l6) + rot(l8) ) 
     &           + c110*( rot(l1) + rot(l3) + rot(l5) + rot(l7) )
     &           + c100*rot(l4) 

             !plan k
             l0 = l0 + inck
             l1 = l1 + inck
             l2 = l2 + inck
             l3 = l3 + inck
             l4 = l4 + inck
             l5 = l5 + inck
             l6 = l6 + inck
             l7 = l7 + inck
             l8 = l8 + inck

             fv =  c110*( rot(l0) + rot(l2) + rot(l6) + rot(l8) ) 
     &           + c100*( rot(l1) + rot(l3) + rot(l5) + rot(l7) )
     &           + c000*rot(l4) 
     &           + fv

             !plan k+1
             l0 = l0 + inck
             l1 = l1 + inck
             l2 = l2 + inck
             l3 = l3 + inck
             l4 = l4 + inck
             l5 = l5 + inck
             l6 = l6 + inck
             l7 = l7 + inck
             l8 = l8 + inck

             fv =  c111*( rot(l0) + rot(l2) + rot(l6) + rot(l8) ) 
     &           + c110*( rot(l1) + rot(l3) + rot(l5) + rot(l7) )
     &           + c100*rot(l4) 
     &           + fv


