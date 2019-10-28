c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine copy_ventijk( ndo, ithread, 
     &        param_int, param_real,
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
c_I    tijk     : vecteur param_int( IO_THREAD)rmale aux facettes des mailles
c_I    vent     : vitesses d'entrainement aux facettes preced.
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    grad    : increment de la solution
c
c***********************************************************************
      implicit none

      INTEGER_E ndo, ithread, ind_dm(6),ind_loop1(6),param_int(0:*)

      REAL_E venti(*),ventj(*),ventk(*),vent_vertex(*)
      REAL_E param_real(0:*)

C Var loc 
      INTEGER_E ind_loop(6),v1ven,v2ven,v3ven,v2in,v3in,
     & i,j,k,l,lij,lxij,inci,incj,inck,lx1,lx2,lx3,l111,l2,l3

#include "FastS/param_solver.h"
#include "FastS/formule_param.h"
#include "FastS/formule_xyz_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"


      ind_loop = ind_loop1

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

             ventk(l)       = .25*( vent_vertex(l111) 
     &                             +vent_vertex(lx2) 
     &                             +vent_vertex(lx3) 
     &                             +vent_vertex(lx1)  )

             ventk(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in)  )

             ventk(l+v3ven) = .25*( vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in)  )
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

             ventj(l+v3ven) = .25*( vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in)  )
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

             venti(l+v3ven) = .25*( vent_vertex(l111+v3in) 
     &                             +vent_vertex(lx2 +v3in) 
     &                             +vent_vertex(lx3 +v3in) 
     &                             +vent_vertex(lx1 +v3in)  )
            enddo
           enddo
          enddo

      ELSEIF(param_int( ITYPVENT ).eq.1) THEN

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
     &                             +vent_vertex(lx1)  )

             ventk(l+v2ven) = .25*( vent_vertex(l111+v2in) 
     &                             +vent_vertex(lx2 +v2in) 
     &                             +vent_vertex(lx3 +v2in) 
     &                             +vent_vertex(lx1 +v2in)  )
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
            enddo
           enddo
          enddo

      ENDIF

      end
