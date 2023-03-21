c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine copy_ventijk( ndo, ithread, 
     &        param_int, param_real,
     &        x, y, z,
     &        ind_dm, ind_loop1,
     &        venti, ventj, ventk, vent_vertex)

c***********************************************************************
c_U   USER : TERRACOL
c
c     ACT
c_A    Appel du calcul des flux explicites
c
c     VAL
c_V    gaz parfait monoespece
c_V    processeur domaine
c_V    steady/unsteady
c
c     INP
c_I    tijk     : vecteur param_int(IO_THREAD) normales aux facettes des mailles
c_I    vent     : vitesses d entrainement aux facettes preced.
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    grad    : increment de la solution
c
c***********************************************************************
      implicit none

      INTEGER_E ndo,ithread,ind_dm(6),ind_loop1(6),param_int(0:*)

      REAL_E venti(*),ventj(*),ventk(*),vent_vertex(*)
      REAL_E param_real(0:*)
      REAL_E x(*),y(*),z(*)

C Var loc 
      INTEGER_E ind_loop(6),v1ven,v2ven,v3ven,v2in,v3in,
     & i,j,k,l,lij,lxij,inci,incj,inck,lx1,lx2,lx3,l111,l2,l3
      REAL_E tkx10, tkx11, tky10, tky11, tkz10, tkz11
      REAL_E vx10, vx11, vy10, vy11, vz10, vz11
      REAL_E Sx, Sy, Sz, O13
      

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_xyz_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

C 0: 1/4 node velocity, 1: ponderated by triangles
#define FORMULA 0

      ind_loop = ind_loop1
      O13 = 1.D0/3.D0

      if(ind_loop1(1).gt.ind_loop1(2)) return 
      if(ind_loop1(3).gt.ind_loop1(4)) return 
      if(ind_loop1(5).gt.ind_loop1(6)) return

      v1ven =   0
      v2ven =   param_int(NDIMDX_VENT)
      v3ven = 2*param_int(NDIMDX_VENT)
      v2in  =   param_int(NDIMDX_XYZ)
      v3in  = 2*param_int(NDIMDX_XYZ)


      if(ind_loop1(1).eq.1) ind_loop(1)= 1-param_int(NIJK_VENT+3)
      if(ind_loop1(3).eq.1) ind_loop(3)= 1-param_int(NIJK_VENT+3)
      if(ind_loop1(5).eq.1) ind_loop(5)= 1-param_int(NIJK_VENT+4)

      if(ind_loop1(2).eq.ind_dm(2)) 
     &  ind_loop(2)=ind_loop1(2)+param_int(NIJK_VENT+3)
      if(ind_loop1(4).eq.ind_dm(4))
     &  ind_loop(4)=ind_loop1(4)+param_int(NIJK_VENT+3)
      if(ind_loop1(6).eq.ind_dm(6)) 
     &  ind_loop(6)=ind_loop1(6)+param_int(NIJK_VENT+4)


      IF(param_int(ITYPVENT).eq.0) THEN

       do k = ind_loop(5), ind_loop(6)
       do j = ind_loop(3), ind_loop(4)

#include    "FastS/Compute/ALE/loopIvent_begin.for"

            !Face k
            inci = 1
            incj = param_int( NIJK_XYZ )

            lx1 = l111 + incj
            lx2 = l111 + inci
            lx3 = lx1  + inci

            ventk(l)       = .25*(  vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

            ventk(l+v2ven) = .25*(  vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

            ventk(l+v3ven) = .25*(  vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in) )

            !formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))
            tkz10=(x(lx1)-x(l111))*(y(lx3)-y(l111))
     &           -(y(lx1)-y(l111))*(x(lx3)-x(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))
            tkz11=(x(lx3)-x(l111))*(y(lx2)-y(l111))
     &           -(y(lx3)-y(l111))*(x(lx2)-x(l111))

            Sx = 0.5*(tkx11 + tkx10)
            Sy = 0.5*(tky11 + tky10)
            Sz = 0.5*(tkz11 + tkz10)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz10 = O13 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx1+v3in)
     &                  +vent_vertex(lx3+v3in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz11 = O13 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx2+v3in)
     &                  +vent_vertex(lx3+v3in))

C            WRITE(*,*) 'kS', Sz
C            WRITE(*,*) 'kref', ventk(l+v3ven)
            IF (ABS(Sx) .GT. 1.e-13) THEN
                ventk(l) = 0.5*(tkx11*vx11 + tkx10*vx10) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN
                ventk(l+v2ven) = 0.5*(tky11*vy11 + tky10*vy10) /Sy
            ENDIF
            IF (ABS(Sz) .GT. 1.e-13) THEN
                ventk(l+v3ven) = 0.5*(tkz11*vz11 + tkz10*vz10) /Sz
            ENDIF
C            WRITE(*,*) 'kmod', ventk(l+v3ven)
#endif
            !Face j
            inci = 1
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + inci
            lx3 = lx1  + inci

            ventj(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

            ventj(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

            ventj(l+v3ven) = .25*( vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in) )

            !formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))
            tkz10=(x(lx1)-x(l111))*(y(lx3)-y(l111))
     &           -(y(lx1)-y(l111))*(x(lx3)-x(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))
            tkz11=(x(lx3)-x(l111))*(y(lx2)-y(l111))
     &           -(y(lx3)-y(l111))*(x(lx2)-x(l111))

            Sx = 0.5*(tkx11 + tkx10)
            Sy = 0.5*(tky11 + tky10)
            Sz = 0.5*(tkz11 + tkz10)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz10 = 013 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx1+v3in)
     &                  +vent_vertex(lx3+v3in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz11 = O13 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx2+v3in)
     &                  +vent_vertex(lx3+v3in))
            
C            WRITE(*,*) 'jS', Sx, Sy
C            WRITE(*,*) 'jref', ventj(l), ventj(l+v2ven)
            IF (ABS(Sx) .GT. 1.e-13) THEN      
                ventj(l) = 0.5*(tkx11*vx11 + tkx10*vx10) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN      
                ventj(l+v2ven) = 0.5*(tky11*vy11 + tky10*vy10) /Sy
            ENDIF
            IF (ABS(Sz) .GT. 1.e-13) THEN      
                ventj(l+v3ven) = 0.5*(tkz11*vz11 + tkz10*vz10) /Sz
            ENDIF
C            WRITE(*,*) 'jmod', ventj(l), ventj(l+v2ven)
#endif

            !Face i
            incj = param_int( NIJK_XYZ )
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + incj
            lx3 = lx1  + incj

            venti(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

            venti(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

            venti(l+v3ven) = .25*( vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in) )

            !formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))
            tkz10=(x(lx1)-x(l111))*(y(lx3)-y(l111))
     &           -(y(lx1)-y(l111))*(x(lx3)-x(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))
            tkz11=(x(lx3)-x(l111))*(y(lx2)-y(l111))
     &           -(y(lx3)-y(l111))*(x(lx2)-x(l111))

            Sx = 0.5*(tkx11 + tkx10)
            Sy = 0.5*(tky11 + tky10)
            Sz = 0.5*(tkz11 + tkz10)            

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz10 = O13 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx1+v3in)
     &                  +vent_vertex(lx3+v3in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))
            vz11 = O13 *(vent_vertex(l111+v3in)
     &                  +vent_vertex(lx2+v3in)
     &                  +vent_vertex(lx3+v3in))
            
C            WRITE(*,*) 'iS', Sx, Sy
C            WRITE(*,*) 'iref', venti(l), venti(l+v2ven)
            IF (ABS(Sx) .GT. 1.e-13) THEN      
                venti(l) = 0.5*(tkx11*vx11 + tkx10*vx10) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN      
                venti(l+v2ven) = 0.5*(tky11*vy11 + tky10*vy10) /Sy
            ENDIF
            IF (ABS(Sz) .GT. 1.e-13) THEN      
                venti(l+v3ven) = 0.5*(tkz11*vz11 + tkz10*vz10) /Sz
            ENDIF
C            WRITE(*,*) 'imod', venti(l), venti(l+v2ven)
#endif

            enddo
           enddo
          enddo

      ELSEIF(param_int(ITYPVENT).eq.1) THEN

       do k = ind_loop(5), ind_loop(6)
       do j = ind_loop(3), ind_loop(4)

#include    "FastS/Compute/ALE/loopIvent_begin.for"

            !Face k
            inci = 1
            incj = param_int( NIJK_XYZ )

            lx1 = l111 + incj
            lx2 = l111 + inci
            lx3 = lx1  + inci

             ventk(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

             ventk(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

            ! formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx3)-z(l111))

            Sx = 0.5*(tkx10 + tkx11)
            Sy = 0.5*(tky10 + tky11)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))

            IF (ABS(Sx) .GT. 1.e-13) THEN            
                ventk(l) = 0.5*(tkx10*vx10 + tkx11*vx11) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN
                ventk(l+v2ven) = 0.5*(tky10*vy10 + tky11*vy11) /Sy
            ENDIF
#endif

            !Face j
            inci = 1
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + inci
            lx3 = lx1  + inci

             ventj(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

             ventj(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

    !formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))

            Sx = 0.5*(tkx10 + tkx11)
            Sy = 0.5*(tky10 + tky11)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))
            IF (ABS(Sx) .GT. 1.e-13) THEN
                ventj(l) = 0.5*(tkx10*vx10 + tkx11*vx11) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN
                ventj(l+v2ven) = 0.5*(tky10*vy10 + tky11*vy11) /Sy
            ENDIF
#endif

            !Face i
            incj = param_int( NIJK_XYZ )
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + incj
            lx3 = lx1  + incj

             venti(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1) )

             venti(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in) )

            !formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))

            Sx = 0.5*(tkx10 + tkx11)
            Sy = 0.5*(tky10 + tky11)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1+v2in)
     &                  +vent_vertex(lx3+v2in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2+v2in)
     &                  +vent_vertex(lx3+v2in))

            IF (ABS(Sx) .GT. 1.e-13) THEN                  
                venti(l) = 0.5*(tkx10*vx10 + tkx11*vx11) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN
                venti(l+v2ven) = 0.5*(tky10*vy10 + tky11*vy11) /Sy
            ENDIF
#endif
            enddo
           enddo
          enddo


      ELSEIF(param_int( ITYPVENT ).eq.3) THEN !2d

       do k = ind_loop(5), ind_loop(6)
       do j = ind_loop(3), ind_loop(4)

#include    "FastS/Compute/ALE/loopIvent_begin.for"

            !Face j
            inci = 1
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + inci
            lx3 = lx1  + inci

             ventj(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1)  )

             ventj(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in)  )

            ! formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))

            Sx = 0.5*(tkx10 + tkx11)
            Sy = 0.5*(tky10 + tky11)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1 +v2in)
     &                  +vent_vertex(lx3 +v2in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2 +v2in)
     &                  +vent_vertex(lx3 +v2in))

            IF (ABS(Sx) .GT. 1.e-13) THEN            
                ventj(l) = 0.5*(tkx10*vx10 + tkx11*vx11) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN
                ventj(l+v2ven) = 0.5*(tky10*vy10 + tky11*vy11) /Sy
            ENDIF
#endif

            !Face i
            incj = param_int( NIJK_XYZ )
            inck = param_int( NIJK_XYZ )*param_int( NIJK_XYZ+1 )

            lx1 = l111 + inck
            lx2 = l111 + incj
            lx3 = lx1  + incj

             venti(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1)  )

             venti(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in)  )

            ! formule de pechier
#if FORMULA == 1
            tkx10=(y(lx1)-y(l111))*(z(lx3)-z(l111))
     &           -(z(lx1)-z(l111))*(y(lx3)-y(l111))
            tky10=(z(lx1)-z(l111))*(x(lx3)-x(l111))
     &           -(x(lx1)-x(l111))*(z(lx3)-z(l111))

            tkx11=(y(lx3)-y(l111))*(z(lx2)-z(l111))
     &           -(z(lx3)-z(l111))*(y(lx2)-y(l111))
            tky11=(z(lx3)-z(l111))*(x(lx2)-x(l111))
     &           -(x(lx3)-x(l111))*(z(lx2)-z(l111))

            Sx = 0.5*(tkx10 + tkx11)
            Sy = 0.5*(tky10 + tky11)

            vx10 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx1)
     &                  +vent_vertex(lx3))
            vy10 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx1 +v2in)
     &                  +vent_vertex(lx3 +v2in))
            
            vx11 = O13 *(vent_vertex(l111)
     &                  +vent_vertex(lx2)
     &                  +vent_vertex(lx3))
            vy11 = O13 *(vent_vertex(l111+v2in)
     &                  +vent_vertex(lx2 +v2in)
     &                  +vent_vertex(lx3 +v2in))
            
            IF (ABS(Sx) .GT. 1.e-13) THEN      
                venti(l) = 0.5*(tkx10*vx10 + tkx11*vx11) /Sx
            ENDIF
            IF (ABS(Sy) .GT. 1.e-13) THEN      
                venti(l+v2ven) = 0.5*(tky10*vy10 + tky11*vy11) /Sy
            ENDIF
#endif
            enddo
           enddo
          enddo

      ENDIF

      end
