c***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cptijk_ale(nd, param_int, ind_loop,
     &                      rot, ti0,tj0,tk0,
     &                      ti,tj,tk)
c***********************************************************************
c_U   USER : ALFEREZ
c
c     ACT
c_A    Calcul des vitesses aux centres des faces, calcul des normales 
c      instationnaires-rotation solide rigide.
c
c     VAL
c_V    Cell center
c
c     INP

c     OUT
c_O    ventk    
c***********************************************************************
      implicit none

#include "FastC/param_solver.h"

      INTEGER_E nd,ind_loop(6), param_int(0:*)

      REAL_E  ti( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tj( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tk( param_int(NDIMDX_MTR) * param_int(NEQ_K ) )
      REAL_E  ti0( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tj0( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) ),
     &        tk0( param_int(NDIMDX_MTR) * param_int(NEQ_K ) )

      REAL_E rot(4,3)
    
c Var Local
      INTEGER_E i,j,k,l1,l2,l3,v1mtr,v2mtr,v3mtr,li,lj,lk
      
#include "FastC/formule_mtr_param.h"

      li   = 0
      lj   = 0
      lk   = 0

      v1mtr =   0
      v2mtr =   param_int(NDIMDX_MTR)
      v3mtr = 2*param_int(NDIMDX_MTR)


      !corection sous domaine pour passer list gradient a liste interface
        if(ind_loop(2).ge.param_int(IJKV)   +1) li   = 1
        if(ind_loop(4).ge.param_int(IJKV+1) +1) lj   = 1
        if(ind_loop(6).ge.param_int(IJKV+2) +1) lk   = 1

      IF(param_int(ITYPZONE).eq.0)THEN !Domaine 3D curvi

      do  k=ind_loop(5),ind_loop(6) 
        do  j=ind_loop(3),ind_loop(4)
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
          do  i=ind_loop(1),ind_loop(2)+li

          l1 = indmtr(i  ,j  ,k ) ! (i  , j  , k  )
          l2 = l1 +  v2mtr
          l3 = l1 +  v3mtr

          ti(l1) =rot(1,1)*ti0(l1) +rot(1,2)*ti0(l2) +rot(1,3)*ti0(l3)
          ti(l2) =rot(2,1)*ti0(l1) +rot(2,2)*ti0(l2) +rot(2,3)*ti0(l3)
          ti(l3) =rot(3,1)*ti0(l1) +rot(3,2)*ti0(l2) +rot(3,3)*ti0(l3)
        enddo
       enddo
      enddo
      do  k=ind_loop(5),ind_loop(6) 
        do  j=ind_loop(3),ind_loop(4)+lj
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
          do  i=ind_loop(1),ind_loop(2)

          l1 = indmtr(i  ,j  ,k     ) ! (i  , j  , k  )
          l2 = l1 +  v2mtr
          l3 = l1 +  v3mtr

          tj(l1) =rot(1,1)*tj0(l1) +rot(1,2)*tj0(l2) +rot(1,3)*tj0(l3)
          tj(l2) =rot(2,1)*tj0(l1) +rot(2,2)*tj0(l2) +rot(2,3)*tj0(l3)
          tj(l3) =rot(3,1)*tj0(l1) +rot(3,2)*tj0(l2) +rot(3,3)*tj0(l3)

        enddo
       enddo
      enddo
      do  k=ind_loop(5),ind_loop(6)+lk
        do  j=ind_loop(3),ind_loop(4)
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
          do  i=ind_loop(1),ind_loop(2)

          l1 = indmtr(i  ,j  ,k     ) ! (i  , j  , k  )
          l2 = l1 +  v2mtr
          l3 = l1 +  v3mtr

          tk(l1) =rot(1,1)*tk0(l1) +rot(1,2)*tk0(l2) +rot(1,3)*tk0(l3)
          tk(l2) =rot(2,1)*tk0(l1) +rot(2,2)*tk0(l2) +rot(2,3)*tk0(l3)
          tk(l3) =rot(3,1)*tk0(l1) +rot(3,2)*tk0(l2) +rot(3,3)*tk0(l3)
        enddo
       enddo
      enddo

      ELSEIF(param_int(ITYPZONE).eq.1.or.param_int(ITYPZONE).eq.3)THEN !Domaine 3D avec une direction homogene k ou 2D:

        !
        !on calcul les metriques du premier plan k
        !
        !
        k=ind_loop(5)
        do  j=ind_loop(3),ind_loop(4)
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
          do  i=ind_loop(1),ind_loop(2)+li

          l1 = indmtr(i  ,j  ,k     ) ! (i  , j  , k  )
          l2 = l1 +  v2mtr
          l3 = l1 +  v3mtr

          ti(l1) =rot(1,1)*ti0(l1) +rot(1,2)*ti0(l2) 
          ti(l2) =rot(2,1)*ti0(l1) +rot(2,2)*ti0(l2) 
        enddo
       enddo
        do  j=ind_loop(3),ind_loop(4)+lj
#ifdef _OPENMP4
CCCC!$OMP simd aligned(ti,tj,tk,ti0,tj0,tk0: CACHELINE)
!$OMP simd 
#else
!DIR$ IVDEP
#endif
          do  i=ind_loop(1),ind_loop(2)

          l1 = indmtr(i  ,j  ,k     ) ! (i  , j  , k  )
          l2 = l1 +  v2mtr
          l3 = l1 +  v3mtr

          tj(l1) =rot(1,1)*tj0(l1) +rot(1,2)*tj0(l2) 
          tj(l2) =rot(2,1)*tj0(l1) +rot(2,2)*tj0(l2) 

        enddo
       enddo

      ENDIF
 
      end
