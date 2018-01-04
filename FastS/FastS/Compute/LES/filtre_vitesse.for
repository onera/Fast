             !plan k-1
             l0  = l  -inci -incj - inck  + var    !(i-1,j-1,k-1)
             l1  = l0 +inci 
             l2  = l1 +inci

             l3  = l0       +incj 
             l4  = l3 +inci 
             l5  = l4 +inci 

             l6  = l3       +incj
             l7  = l4       +incj 
             l8  = l5       +incj 
             
             fv =  c111*( rop(l0) + rop(l2) + rop(l6) + rop(l8) ) 
     &           + c110*( rop(l1) + rop(l3) + rop(l5) + rop(l7) )
     &           + c100*rop(l4) 

c             if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,4f19.15,3i4)')'k0i1',
c     &     rop(l0) , rop(l2) , rop(l6) , rop(l8)        ,l-lij,j,k
c        write(*,'(a,4f19.15)')'k0i2',
c     &     rop(l1) , rop(l3) , rop(l5) , rop(l7)        
c        write(*,'(a,2f19.15)')'k0i3',  rop(l4) ,fv
c             endif

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

             fv =  c110*( rop(l0) + rop(l2) + rop(l6) + rop(l8) ) 
     &           + c100*( rop(l1) + rop(l3) + rop(l5) + rop(l7) )
     &           + c000*rop(l4) 
     &           + fv

c            if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,4f19.15,3i4)')'k1i1',
c     &     rop(l0) , rop(l2) , rop(l6) , rop(l8)        ,l-lij,j,k
c        write(*,'(a,4f19.15)')'k1i2',
c     &     rop(l1) , rop(l3) , rop(l5) , rop(l7)        
c        write(*,'(a,2f19.15)')'k1i3',  rop(l4) ,fv
c             endif

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

             fv =  c111*( rop(l0) + rop(l2) + rop(l6) + rop(l8) ) 
     &           + c110*( rop(l1) + rop(l3) + rop(l5) + rop(l7) )
     &           + c100*rop(l4) 
     &           + fv

c            if(k.eq.20.and.j.le.2.and.l-lij.le.3) then
c        write(*,'(a,4f19.15,3i4)')'k2i1',
c     &     rop(l0) , rop(l2) , rop(l6) , rop(l8)        ,l-lij,j,k
c        write(*,'(a,4f19.15)')'k2i2',
c     &     rop(l1) , rop(l3) , rop(l5) , rop(l7)        
c        write(*,'(a,2f19.15)')'k2i3',  rop(l4) ,fv
c             endif
