c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 35 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bvbs_inflow_supersonic(idir,lrhs, neq_mtr,
     &                                  param_int, ind_loop,
     &                                  param_real, c4,c5,c6,
     &                                 ventijk, tijk, rop, state)
c_U   USER : PECHIER
c
c     ACT
c     Paroi glissante
c     VAL
c_V    Optimisation NEC
c
c     COM
c_C    MODIFIER bvas3.f (passage de ird1,ird2a,ird2b ?????)
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E idir,lrhs, neq_mtr, ind_loop(6), param_int(0:*)

      REAL_E rop    (param_int(NDIMDX     ), param_int(NEQ)      )
      REAL_E ventijk(param_int(NDIMDX_VENT), param_int(NEQ_VENT) )
      REAL_E tijk   (param_int(NDIMDX_MTR ), neq_mtr             )
      REAL_E state(param_int(NEQ))
      REAL_E c4,c5,c6, param_real(0:*)

C Var local
      INTEGER_E  l,l0,incj,inck,lij,lr,lp,i,j,k,l1,l2,inci
      REAL_E ro,u,v,w,t,nut,c0,c1,c2,c3

#include "FastS/formule_param.h"


      !!! premiere rangee = state
      !!! rangee suivante pour que : Phi_L + Phi_R = 2 * state

      ro    = state(1)
      u     = state(2)/state(1)
      v     = state(3)/state(1)
      w     = state(4)/state(1)
      t     = ( state(5) -0.5*ro*(u*u+v*v+w*w) ) /(ro*param_real(CVINF))

      if (param_int(NEQ).eq.6) nut = state(6)/state(1)

      !write(*,'(a,5f18.12)')'state', ro,u,v,w,t

      c0   = 1./c6
      c1   =-(c4 + c5)*c0
      c2   =- c6*c0
      c3   = (2.- c5- c4)*c0


      IF (idir.eq.1) THEN

          if(param_int(NEQ).eq.5) then

             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

               l    = inddm(  ind_loop(2) , j,  k ) 
#include       "FastS/BC/BCInflowSupersonic_firstrank.for"

               l0   = l
               l1   = l + 1
               l2   = l + 2
               do i = ind_loop(1), ind_loop(2)-1

                   l    = inddm( i , j,  k ) 
#include           "FastS/BC/BC_nextrank.for"
               enddo
             enddo
             enddo

          else

             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)

               l    = inddm(  ind_loop(2) , j,  k ) 
#include       "FastS/BC/BCInflowSupersonic_firstrank_SA.for"

               l0   = l
               l1   = l + 1
               l2   = l + 2
               do i = ind_loop(1), ind_loop(2)-1

                   l    = inddm( i , j,  k ) 
#include           "FastS/BC/BC_nextrank_SA.for"
               enddo
             enddo
             enddo

          endif !param_int(NEQ)

      ELSEIF (idir.eq.2) THEN

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(1)    , j,  k ) 
#include      "FastS/BC/BCInflowSupersonic_firstrank.for"

               l0   = l
               l1   = l - 1
               l2   = l - 2
               do i = ind_loop(1)+1, ind_loop(2)

                   l    = inddm( i , j,  k ) 
#include           "FastS/BC/BC_nextrank.for"
               enddo
             enddo
             enddo

          else

             do  k = ind_loop(5), ind_loop(6)
             do  j = ind_loop(3), ind_loop(4)

               l    = inddm( ind_loop(1)    , j,  k ) 
#include      "FastS/BC/BCInflowSupersonic_firstrank_SA.for"

               l0   = l
               l1   = l - 1
               l2   = l - 2
               do i = ind_loop(1)+1, ind_loop(2)

                   l    = inddm( i , j,  k ) 
#include           "FastS/BC/BC_nextrank_SA.for"
               enddo
             enddo
             enddo

          endif !param_int(NEQ)


      ELSEIF (idir.eq.3) THEN

          incj = param_int(NIJK)

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)

                  lij =  inddm( ind_loop(1),  ind_loop(4)    , k )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank.for"
                  enddo !i

                do  j = ind_loop(3), ind_loop(4)-1

                  lij =       inddm( ind_loop(1),     j        , k )
                  lr  = lij - inddm( ind_loop(1),  ind_loop(4) , k )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                     l0   = l -lr
                     l1   = l0 +   incj
                     l2   = l0 + 2*incj
#include             "FastS/BC/BC_nextrank.for"
                  enddo!i
                enddo !j
             enddo !k

          else

             do  k = ind_loop(5), ind_loop(6)

                  lij =  inddm( ind_loop(1),  ind_loop(4)    , k )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank_SA.for"
                  enddo !i

                  do  j = ind_loop(3), ind_loop(4)-1

                    lij =       inddm( ind_loop(1),     j        , k )
                    lr  = lij - inddm( ind_loop(1),  ind_loop(4) , k )
!DEC$ IVDEP
                    do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                       l0   = l -lr
                       l1   = l0 +   incj
                       l2   = l0 + 2*incj
#include               "FastS/BC/BC_nextrank_SA.for"
                    enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)

      ELSEIF (idir.eq.4) THEN

          incj = -param_int(NIJK)

          if(param_int(NEQ).eq.5) then

             do  k = ind_loop(5), ind_loop(6)

                  lij =  inddm( ind_loop(1),  ind_loop(3)    , k )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank.for"
                  enddo !i

                  do  j = ind_loop(3)+1, ind_loop(4)

                    lij =       inddm( ind_loop(1),     j        , k )
                    lr  = lij - inddm( ind_loop(1),  ind_loop(3) , k )
!DEC$ IVDEP
                    do l = lij, lij + ind_loop(2) - ind_loop(1)
     
                       l0   = l -lr
                       l1   = l0 +   incj
                       l2   = l0 + 2*incj
#include               "FastS/BC/BC_nextrank.for"
                    enddo!i
                enddo !j
             enddo !k

          else

             do  k = ind_loop(5), ind_loop(6)

                  lij =  inddm( ind_loop(1),  ind_loop(3)    , k )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank_SA.for"
                  enddo !i

                  do  j = ind_loop(3)+1, ind_loop(4)

                    lij =       inddm( ind_loop(1),     j        , k )
                    lr  = lij - inddm( ind_loop(1),  ind_loop(3) , k )
!DEC$ IVDEP
                    do l = lij, lij + ind_loop(2) - ind_loop(1)
     
                       l0   = l -lr
                       l1   = l0 +   incj
                       l2   = l0 + 2*incj
#include               "FastS/BC/BC_nextrank_SA.for"
                    enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)


      ELSEIF (idir.eq.5) THEN


          inck = param_int(NIJK)*param_int(NIJK+1)

          if(param_int(NEQ).eq.5) then

             do  j = ind_loop(3), ind_loop(4)

                  lij =  inddm( ind_loop(1),  ind_loop(4), ind_loop(6) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank.for"
                  enddo !i
             enddo !j

             do  k = ind_loop(5), ind_loop(6)-1

                do  j = ind_loop(3), ind_loop(4)

                  lij =       inddm( ind_loop(1),  j  ,    k        )
                  lr  = lij - inddm( ind_loop(1),  j  , ind_loop(6) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                     l0   = l -lr
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
                  enddo!i
                enddo !j
             enddo !k

          else

             do  j = ind_loop(3), ind_loop(4)

                  lij =  inddm( ind_loop(1),  ind_loop(4), ind_loop(6) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank_SA.for"
                  enddo !i
             enddo !j

             do  k = ind_loop(5), ind_loop(6)-1

                do  j = ind_loop(3), ind_loop(4)

                  lij =       inddm( ind_loop(1),  j  ,    k        )
                  lr  = lij - inddm( ind_loop(1),  j  , ind_loop(6) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                     l0   = l -lr
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank_SA.for"
                  enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)

      ELSE 

          inck = -param_int(NIJK)*param_int(NIJK+1)

          if(param_int(NEQ).eq.5) then

             do  j = ind_loop(3), ind_loop(4)

                  lij =  inddm( ind_loop(1),  ind_loop(4), ind_loop(5) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank.for"
                  enddo !i
             enddo !j

             do  k = ind_loop(5)+1, ind_loop(6)

                do  j = ind_loop(3), ind_loop(4)

                  lij =       inddm( ind_loop(1),  j  ,    k        )
                  lr  = lij - inddm( ind_loop(1),  j  , ind_loop(5) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                     l0   = l -lr
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank.for"
                  enddo!i
                enddo !j
             enddo !k

          else

             do  j = ind_loop(3), ind_loop(4)

                  lij =  inddm( ind_loop(1),  ind_loop(4), ind_loop(5) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
#include           "FastS/BC/BCInflowSupersonic_firstrank_SA.for"
                  enddo !i
             enddo !j

             do  k = ind_loop(5)+1, ind_loop(6)

                do  j = ind_loop(3), ind_loop(4)

                  lij =       inddm( ind_loop(1),  j  ,    k        )
                  lr  = lij - inddm( ind_loop(1),  j  , ind_loop(5) )
!DEC$ IVDEP
                  do l = lij, lij + ind_loop(2) - ind_loop(1)
   
                     l0   = l -lr
                     l1   = l0 +   inck
                     l2   = l0 + 2*inck
#include             "FastS/BC/BC_nextrank_SA.for"
                  enddo!i
                enddo !j
             enddo !k

          endif !param_int(NEQ)

      ENDIF !idir


      END
