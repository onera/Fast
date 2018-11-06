c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 39 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bflinj1(ndom, idir, mobile_coef, neq_mtr,
     &                   param_int, param_real, incijk, ind_loop,
     &                   rop, drodm, tijk, ventijk)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  bflco0k.sf/pechier1/e1/i1
c
c_DC  DATE_C : Octobre 1990 -- GUILLEN / DORMIEUX
c
c     HISTORIQUE
c_H    Iteration     1 ___ 1997-06-05 15:34:15  pechier
c_H    Issu de bflco0k.sf/lepape1/e1/i4
c
c_U   USER : GUILLEN
c
c     ACT
c_A    Calcul de flux frontieres aux interfaces solides du domaine,
c_A    la pression est obtenue par relations de compatibilite.
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    processeur domaine
c
c     INP
c_I    ndom   : numero du domaine
c_I    idir   : direction du traitement
c_I    ndimdx : dimension 3D
c_I    tijk   : tableaux des normales
c_I    flu    : tableau des flux explicites
c_I    q      : variables d'etats aux interfaces
cbc_type.eq.3
c     OUT
c
c     I/O
c_/    flu    : increment des variables conservatives
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,neq_mtr,incijk, idir, ind_loop(6), param_int(0:*)

      REAL_E   rop( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E drodm( param_int(NDIMDX)     * param_int(NEQ)     )
      REAL_E ventijk( param_int(NDIMDX_VENT)* param_int(NEQ_VENT))

      REAL_E  tijk( param_int(NDIMDX_MTR) * param_int(NEQ_IJ) )

      REAL_E param_real(0:*), mobile_coef

C var loc
      INTEGER_E im,jm,km,ijkm,l,l0,iadrf,i,j,k,lj,ic,jc,kc,
     & kc_vent,v1,v2,v3,v4,v5,v6,vmtr,vven,shift,lt,lven,lij,ltij,lvij

      REAL_E p,r,u,v,w,qen,ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale,
     & ck_mtr_vent, u_int,sens,rgp,flu1,flu2,flu3,flu4,flu5,flu6,
     & tcx,tcy,tcz,gam2, h
 
#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"
#include "FastS/formule_vent_param.h"

      rgp  = param_real( CVINF )*(param_real( GAMMA )-1.)  !Cv(gama-1)= R (gas parfait)
      gam2 = param_real( GAMMA )/ (param_real( GAMMA )-1.)

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)

      vmtr  =   param_int(NDIMDX_MTR)
      vven  =   param_int(NDIMDX_VENT)

      sens  = 1.
      shift = 0
      if(mod(idir,2).eq.0) then
        sens  =-1.
        shift = 1
      endif

      shift =incijk*shift

      !!  a mettre dans Fast
      call shape_tab_mtr(neq_mtr, param_int, idir,
     &                   ic,jc,kc,kc_vent,
     &                   ci_mtr,cj_mtr,ck_mtr,ck_vent,c_ale)

      c_ale = c_ale*mobile_coef


      !correction monoindice
      ic      = ic -1
      jc      = jc -1
      kc      = kc -1
      kc_vent = kc_vent -1

      IF(param_int(NEQ).eq.5) THEN

         DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)
         DO i = ind_loop(1), ind_loop(2)

          l     = inddm(  i, j, k)
          lt    = indmtr( i, j, k)
          lven  = indven( i, j, k)

          iadrf = l  - incijk
          l0    = l  - shift

          tcx = tijk(lt +vmtr*ic)*ci_mtr
          tcy = tijk(lt +vmtr*jc)*cj_mtr
          tcz = tijk(lt +vmtr*kc)*ck_mtr

          qen =(  ventijk(lven              )*tcx
     &           +ventijk(lven +vven        )*tcy
     &           +ventijk(lven +vven*kc_vent)*tcz*ck_vent
     &         )*c_ale

          r     = 0.5*(rop(l+v1)+rop(iadrf+v1))
          u     = 0.5*(rop(l+v2)+rop(iadrf+v2))
          v     = 0.5*(rop(l+v3)+rop(iadrf+v3))
          w     = 0.5*(rop(l+v4)+rop(iadrf+v4))

          !determination vitesse normale interface
          u_int= tcx*u +tcy*v +tcz*w -qen


          r     =rop(iadrf+v1)
          u     =rop(iadrf+v2)
          v     =rop(iadrf+v3)
          w     =rop(iadrf+v4)

          p     = rop(iadrf+v5)*rop(iadrf+v1)*rgp

          h     = gam2*p + r*(u*u+v*v+w*w)*0.5

          p=0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

          flu1=            u_int*r
          flu2= tcx * p  + u_int*r*u
          flu3= tcy * p  + u_int*r*v
          flu4= tcz * p  + u_int*r*w
          flu5=            u_int*h    + p*qen
#include  "FastS/Compute/assemble_drodm_corr.for"
       enddo
       enddo
       enddo

      ELSE

         !write(*,*)'loop',ind_loop(1),incijk,shift
         DO k = ind_loop(5), ind_loop(6)
         DO j = ind_loop(3), ind_loop(4)
         DO i = ind_loop(1), ind_loop(2)

          l     = inddm(  i, j, k)
          lt    = indmtr( i, j, k)
          lven  = indven( i, j, k)

          iadrf = l  - incijk
          l0    = l  - shift

          tcx = tijk(lt +vmtr*ic)*ci_mtr
          tcy = tijk(lt +vmtr*jc)*cj_mtr
          tcz = tijk(lt +vmtr*kc)*ck_mtr

          qen =(  ventijk(lven              )*tcx
     &           +ventijk(lven +vven        )*tcy
     &           +ventijk(lven +vven*kc_vent)*tcz*ck_vent
     &         )*c_ale

          r     = 0.5*(rop(l+v1)+rop(iadrf+v1))
          u     = 0.5*(rop(l+v2)+rop(iadrf+v2))
          v     = 0.5*(rop(l+v3)+rop(iadrf+v3))
          w     = 0.5*(rop(l+v4)+rop(iadrf+v4))

          !determination vitesse normale interface
          u_int= tcx*u +tcy*v +tcz*w -qen

          r     =rop(iadrf+v1)
          u     =rop(iadrf+v2)
          v     =rop(iadrf+v3)
          w     =rop(iadrf+v4)

          p     = rop(iadrf+v5)*rop(iadrf+v1)*rgp

          h     = gam2*p + r*(u*u+v*v+w*w)*0.5

          p=0.5*(rop(l+v5)*rop(l+v1)+rop(iadrf+v5)*rop(iadrf+v1))*rgp

          flu1=            u_int*r
          flu2= tcx * p  + u_int*r*u
          flu3= tcy * p  + u_int*r*v
          flu4= tcz * p  + u_int*r*w
          flu5=            u_int*h    + p*qen
          flu6=            u_int*r*rop(iadrf+v6)
#include    "FastS/Compute/SA/assemble_drodm_corr.for"
       enddo
       enddo
       enddo
      ENDIF

      end
