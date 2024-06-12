c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine move_to_temp( nd, ndom, param_int, rop, temp)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  mjdro3.sf/pechier1/e1/i1
c
c     ACT
c_A    Mise a jour de drodm apres calcul des flux.
c
c     INP
c_I    ndom   : numero-utilisateur du domaine
c_I    ipara  : direction du calcul des flux
c_I    neq    : nombre d equations du systeme
c_I    ndimdx : nombre de points maximal dans un domaine
c_I    flu    : flux aux interfaces
c
c     OUT
c
c     I/O
c_/    drodm  : increment des variables conservatives
c***********************************************************************
      implicit none

#include "Fast/param_solver.h"

      INTEGER_E nd, ndom, param_int(0:*)

      REAL_E rop(param_int(NDIMDX), param_int(NEQ)  )
      REAL_E temp(param_int(NDIMDX), param_int(NEQ) )

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,b,ind,ni,ninj,p,l2,
     & cb( param_int(NEQ_LBM), 3)
      INTEGER_E l1,l3,l4,l5,l6, l_alt, i_idx, l7,l8,l9,l10,l11,l12
      INTEGER_E im2,ip2,im1,ip1,jm2,jp2,jm1,jp1,km2,km1,kp2,kp1
      INTEGER_E ng, nx, ny, nz, ot

C    adresse point courant pour tableau de la taille d'un domaine
       INTEGER_E inddm, i_1,j_1,k_1

       inddm(i_1,j_1,k_1) = 1
     &     + (i_1+param_int(NIJK+3)-1)
     &     + (j_1+param_int(NIJK+3)-1)*param_int(NIJK)
     &     + (k_1+param_int(NIJK+4)-1)*param_int(NIJK)*param_int(NIJK+1)

      ni  = param_int(NIJK)
      ninj= param_int(NIJK)*param_int(NIJK+1)

      nx = param_int(IJKV  )
      ny = param_int(IJKV+1)
      nz = param_int(IJKV+2)
      ng = param_int(NIJK+3)

      print*, 'Bonjour je passe dans move'
      print*, (nz+2*ng)*(ny+2*ng)*(nx+2*ng)

           do  k = 1+ng, nz+ng
            do  j = 1+ng, ny+ng
              do i = 1+ng, nx+ng

                l = inddm(i,j,k)
                l1 = l + 1
                l2 = l - 1
                l3 = l + ni
                l4 = l - ni
                l5 = l + ninj
                l6 = l - ninj
                l7 = l + 2
                l8 = l - 2
                l9 = l + 2*ni
                l10 = l - 2*ni
                l11 = l + 2*ninj
                l12 = l - 2*ninj

                do ot = 1, param_int(NEQ)

                    temp(l, ot) = rop(l ,ot)
                    temp(l1,ot) = rop(l1,ot)
                    temp(l2,ot) = rop(l2,ot)
                    temp(l3,ot) = rop(l3,ot)
                    temp(l4,ot) = rop(l4,ot)
                    temp(l5,ot) = rop(l5,ot)
                    temp(l6,ot) = rop(l6,ot)
                    temp(l7,ot) = rop(l7,ot)
                    temp(l8,ot) = rop(l8,ot)
                    temp(l9,ot) = rop(l9,ot)
                    temp(l10,ot) = rop(l10,ot)
                    temp(l11,ot) = rop(l11,ot)
                    temp(l12,ot) = rop(l12,ot)

                enddo

              enddo
            enddo
           enddo


      end
