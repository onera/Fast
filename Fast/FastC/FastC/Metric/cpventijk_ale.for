***********************************************************************
c     $Date: 2010-12-22 11:30:40 +0100 (Wed, 22 Dec 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine cpventijk_ale(ndom, param_int,param_real,ind_loop_glob,
     &                         x,y,z, rot,
     &                         venti,ventj,ventk)
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

      INTEGER_E ndom,ind_loop_glob(6),param_int(0:*)

      REAL_E x(param_int(NDIMDX_XYZ)),y(param_int(NDIMDX_XYZ)),
     &       z(param_int(NDIMDX_XYZ))

      REAL_E venti(param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E ventj(param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E ventk(param_int(NDIMDX_VENT)* param_int(NEQ_VENT))
      REAL_E rot(4,3), param_real(0:*)
    
c Var Local
      INTEGER_E i,j,k,inci,incj,inck,li,lj,lk,v1ven,v2ven,v3ven,
     & l,l111,l121,l221,l211,l112,l122,l212,l1,l2,l3,lxij,ind_loop(6),
     & lv,lij,ltij,lx,lvij, lt, lvo
      REAL_E cax,cay,caz,cenix,ceniy,ceniz,cenjx,cenjy,cenjz,
     & cenkx,cenky,cenkz,vtrans(3)
      
#include "FastC/formule_param.h"
#include "FastC/formule_xyz_param.h"
#include "FastC/formule_mtr_param.h"
#include "FastC/formule_vent_param.h"

      vtrans= 0.

      li    = 0
      lj    = 0
      lk    = 0
      v1ven = 0
      v2ven =   param_int(NDIMDX_VENT)
      v3ven = 2*param_int(NDIMDX_VENT)

      !corection sous domaine pour passer list gradient a liste interface
      if(ind_loop_glob(2).ge.param_int(IJKV)   +1) li   = 1
      if(ind_loop_glob(4).ge.param_int(IJKV+1) +1) lj   = 1
      if(ind_loop_glob(6).ge.param_int(IJKV+2) +1) lk   = 1

      inci = 1
      incj = param_int(NIJK_XYZ)
      inck = param_int(NIJK_XYZ)*param_int(NIJK_XYZ+1)

      IF(param_int(ITYPVENT).eq.0) THEN !mvt 3dfull

        !!
        !!! Facette I
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(2)= ind_loop(2)+li
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            l111= lx            ! x(i , j  , k  )
            l121= l111 + incj   ! x(i , j+1, k  )
            l112= l111 + inck   ! x(i , j  , k+1)
            l122= l121 + inck   ! x(i , j+1, k+1)

            cenix = .25*( x(l111) + x(l121) + x(l112) + x(l122) )
            ceniy = .25*( y(l111) + y(l121) + y(l112) + y(l122) )
            ceniz = .25*( z(l111) + z(l121) + z(l112) + z(l122) )

            cax =  rot(1,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(ceniz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(ceniz - param_real(ROT_CENTER+2))

            caz =  rot(3,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(3,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(3,3)*(ceniz - param_real(ROT_CENTER+2))

            venti(l1) = vtrans(1) + rot(4,2)*caz - rot(4,3)*cay
            venti(l2) = vtrans(2) + rot(4,3)*cax - rot(4,1)*caz
            venti(l3) = vtrans(3) + rot(4,1)*cay - rot(4,2)*cax
#include "FastC/HPC_LAYER/loop_end.for"
        !!
        !!! Facette J
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(4)= ind_loop(4)+lj
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            l111= lx                        ! x(i  , j  , k  )
            l112= l111 + inck               ! x(i  , j  , k+1)
            l211= l111 + inci               ! x(i+1, j  , k  )
            l212= l112 + inci               ! x(i+1, j  , k+1)

            cenjx = .25*( x(l111) + x(l211) + x(l112) + x(l212) )
            cenjy = .25*( y(l111) + y(l211) + y(l112) + y(l212) )
            cenjz = .25*( z(l111) + z(l211) + z(l112) + z(l212) )

            cax =  rot(1,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenjz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenjz - param_real(ROT_CENTER+2))

            caz =  rot(3,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(3,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(3,3)*(cenjz - param_real(ROT_CENTER+2))

            ventj(l1) = vtrans(1) + rot(4,2)*caz - rot(4,3)*cay
            ventj(l2) = vtrans(2) + rot(4,3)*cax - rot(4,1)*caz
            ventj(l3) = vtrans(3) + rot(4,1)*cay - rot(4,2)*cax
#include "FastC/HPC_LAYER/loop_end.for"
          !!
          !!! Facette K
          !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(6)= ind_loop(6)+lk
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            l111= lx                        ! x(i  , j  , k  )
            l121= l111 + incj               ! x(i  , j+1, k  )
            l211= l111 + inci               ! x(i+1, j  , k  )
            l221= l121 + inci               ! x(i+1, j+1, k  )

            cenkx = .25*( x(l111) + x(l211) + x(l121) + x(l221) )
            cenky = .25*( y(l111) + y(l211) + y(l121) + y(l221) )
            cenkz = .25*( z(l111) + z(l211) + z(l121) + z(l221) )

            cax =  rot(1,1)*(cenkx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(cenky - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenkz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenkx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(cenky - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenkz - param_real(ROT_CENTER+2))

            caz =  rot(3,1)*(cenkx - param_real(ROT_CENTER  ))
     &           + rot(3,2)*(cenky - param_real(ROT_CENTER+1))
     &           + rot(3,3)*(cenkz - param_real(ROT_CENTER+2))

            ventk(l1) = vtrans(1) + rot(4,2)*caz - rot(4,3)*cay
            ventk(l2) = vtrans(2) + rot(4,3)*cax - rot(4,1)*caz
            ventk(l3) = vtrans(3) + rot(4,1)*cay - rot(4,2)*cax
#include "FastC/HPC_LAYER/loop_end.for"

      ELSEIF(param_int(ITYPVENT).eq.1) THEN !mvt 2D (plan (x,y): rot Z pur)

        !!
        !!! Facette I
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(2)= ind_loop(2)+li
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven

            l111= lx                        ! x(i  , j  , k  )
            l121= l111 + incj               ! x(i  , j+1, k  )
            l112= l111 + inck               ! x(i  , j  , k+1)
            l122= l121 + inck               ! x(i  , j+1, k+1)

            cenix = .25*( x(l111) + x(l121) + x(l112) + x(l122) )
            ceniy = .25*( y(l111) + y(l121) + y(l112) + y(l122) )
            ceniz = .25*( z(l111) + z(l121) + z(l112) + z(l122) )

            cax =  rot(1,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(ceniz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(ceniz - param_real(ROT_CENTER+2))

            venti(l1) = vtrans(1) - rot(4,3)*cay
            venti(l2) = vtrans(2) + rot(4,3)*cax
#include "FastC/HPC_LAYER/loop_end.for"
        !!
        !!! Facette J
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(4)= ind_loop(4)+lj
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            l111= indcg(i  ,j  ,k )         ! x(i  , j  , k  )
            l112= l111 + inck               ! x(i  , j  , k+1)
            l211= l111 + inci               ! x(i+1, j  , k  )
            l212= l112 + inci               ! x(i+1, j  , k+1)

            cenjx = .25*( x(l111) + x(l211) + x(l112) + x(l212) )
            cenjy = .25*( y(l111) + y(l211) + y(l112) + y(l212) )
            cenjz = .25*( z(l111) + z(l211) + z(l112) + z(l212) )

            cax =  rot(1,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenjz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenjz - param_real(ROT_CENTER+2))

            ventj(l1) = vtrans(1) - rot(4,3)*cay
            ventj(l2) = vtrans(2) + rot(4,3)*cax
#include "FastC/HPC_LAYER/loop_end.for"
          !!
          !!! Facette K
          !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(6)= ind_loop(6)+lk
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven

            l111= lx                        ! x(i  , j  , k  )
            l121= l111 + incj               ! x(i  , j+1, k  )
            l211= l111 + inci               ! x(i+1, j  , k  )
            l221= l121 + inci               ! x(i+1, j+1, k  )

            cenkx = .25*( x(l111) + x(l211) + x(l121) + x(l221) )
            cenky = .25*( y(l111) + y(l211) + y(l121) + y(l221) )
            cenkz = .25*( z(l111) + z(l211) + z(l121) + z(l221) )

            cax =  rot(1,1)*(cenkx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(cenky - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenkz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenkx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(cenky - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenkz - param_real(ROT_CENTER+2))

            ventk(l1) = vtrans(1) - rot(4,3)*cay
            ventk(l2) = vtrans(2) + rot(4,3)*cax
#include "FastC/HPC_LAYER/loop_end.for"

      ELSEIF(param_int(ITYPVENT).eq.2) THEN !translation

         ind_loop(:)= ind_loop_glob(:)
         ind_loop(2)= ind_loop(2)+li
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            venti(l1) = vtrans(1) 
            venti(l2) = vtrans(2) 
            venti(l3) = vtrans(3) 
#include "FastC/HPC_LAYER/loop_end.for"

      ELSE !  dom 2d


        !!
        !!! Facette I
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(2)= ind_loop(2)+li
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven

            l111= lx                        ! x(i  , j  , k  )
            l121= l111 + incj               ! x(i  , j+1, k  )
            l112= l111 + inck               ! x(i  , j  , k+1)
            l122= l121 + inck               ! x(i  , j+1, k+1)

            cenix = .25*( x(l111) + x(l121) + x(l112) + x(l122) )
            ceniy = .25*( y(l111) + y(l121) + y(l112) + y(l122) )
            ceniz = .25*( z(l111) + z(l121) + z(l112) + z(l122) )

            cax =  rot(1,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(ceniz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenix - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(ceniy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(ceniz - param_real(ROT_CENTER+2))

            venti(l1) = vtrans(1) - rot(4,3)*cay
            venti(l2) = vtrans(2) + rot(4,3)*cax
#include "FastC/HPC_LAYER/loop_end.for"
        !!
        !!! Facette J
        !!
         ind_loop(:)= ind_loop_glob(:)
         ind_loop(4)= ind_loop(4)+lj
#include "FastC/HPC_LAYER/loop_ale_begin.for"
            l1 = lv ! vent(i, j, k, 1)
            l2 = l1 +  v2ven
            l3 = l1 +  v3ven

            l111= lx                        ! x(i  , j  , k  )
            l112= l111 + inck               ! x(i  , j  , k+1)
            l211= l111 + inci               ! x(i+1, j  , k  )
            l212= l112 + inci               ! x(i+1, j  , k+1)

            cenjx = .25*( x(l111) + x(l211) + x(l112) + x(l212) )
            cenjy = .25*( y(l111) + y(l211) + y(l112) + y(l212) )
            cenjz = .25*( z(l111) + z(l211) + z(l112) + z(l212) )

            cax =  rot(1,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(1,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(1,3)*(cenjz - param_real(ROT_CENTER+2))
 
            cay =  rot(2,1)*(cenjx - param_real(ROT_CENTER  ))
     &           + rot(2,2)*(cenjy - param_real(ROT_CENTER+1))
     &           + rot(2,3)*(cenjz - param_real(ROT_CENTER+2))

            ventj(l1) = vtrans(1) - rot(4,3)*cay
            ventj(l2) = vtrans(2) + rot(4,3)*cax
#include "FastC/HPC_LAYER/loop_end.for"

      ENDIF !translation/rotation


      end
