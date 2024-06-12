c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine filtrage_5( nd, ndom, param_int, rop, temp)
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
      INTEGER_E ng, nx, ny, nz

      REAL_E sigma_f, d0, d1, d2
      REAL_E f_ro, f_rou, f_rov, f_row
      REAL_E fili, filj, filk

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

      ! sigma_f = 0.2
      ! d0 =  6./16.
      ! d1 = -4./16.
      ! d2 =  1./16.

      sigma_f = 0.01
      d0 =  1./2.
      d1 = -1./4.
      d2 =  0.

      print*, 'Bonjour je passe dans filtre'
      print*, (nz+2*ng)*(ny+2*ng)*(nx+2*ng)

           do  k = 1+ng, nz+ng
            do  j = 1+ng, ny+ng
              do i = 1+ng, nx+ng

                l = inddm(i,j,k)
                im2 = l-2
                ip2 = l+2
                im1 = l-1
                ip1 = l+1
                jm2 = l-2*ni
                jp2 = l+2*ni
                jm1 = l-ni
                jp1 = l+ni
                km2 = l-2*ninj
                kp2 = l+2*ninj
                km1 = l-ninj
                kp1 = l+ninj

                f_ro = temp(l,1)
                f_rou = temp(l,1)*temp(l,2)
                f_rov = temp(l,1)*temp(l,3)
                f_row = temp(l,1)*temp(l,4)

                fili = 1.
                filj = 1.
                filk = 1.

                !i_idx = mod(l,param_int(NIJK))-param_int(NIJK+3)
                if (i==1+ng) then
                  fili = 0.
                endif
                if (i==nx+ng) then
                  fili = 0.
                endif
                if (j==1+ng) then
                  filj = 0.
                endif
                if (j==ny+ng) then
                  filj = 0.
                endif
                if (k==1+ng) then
                  filk = 0.
                endif
                if (k==nz+ng) then
                  filk = 0.
                endif

                ! selon i
                f_ro =f_ro-fili*sigma_f*(d2*(temp(im2,1)+temp(ip2,1))
     &                           +d1*(temp(im1,1)+temp(ip1,1))
     &                           +d0*temp(l,1))
                f_rou =f_rou-fili*sigma_f*(d2*(temp(im2,1)*temp(im2,2)
     &                                   +temp(ip2,1)*temp(ip2,2))
     &                               +d1*(temp(im1,1)*temp(im1,2)
     &                                   +temp(ip1,1)*temp(ip1,2))
     &                               +d0*temp(l,1)*temp(l,2))
                f_rov =f_rov-fili*sigma_f*(d2*(temp(im2,1)*temp(im2,3)
     &                                     +temp(ip2,1)*temp(ip2,3))
     &                                 +d1*(temp(im1,1)*temp(im1,3)
     &                                     +temp(ip1,1)*temp(ip1,3))
     &                                 +d0*temp(l,1)*temp(l,3))
                f_row =f_row-fili*sigma_f*(d2*(temp(im2,1)*temp(im2,4)
     &                                     +temp(ip2,1)*temp(ip2,4))
     &                                 +d1*(temp(im1,1)*temp(im1,4)
     &                                     +temp(ip1,1)*temp(ip1,4))
     &                                  +d0*temp(l,1)*temp(l,4))

                ! selon j
                f_ro =f_ro-filj*sigma_f*(d2*(temp(jm2,1)+temp(jp2,1))
     &                             +d1*(temp(jm1,1)+temp(jp1,1))
     &                             +d0*temp(l,1))
                f_rou =f_rou-filj*sigma_f*(d2*(temp(jm2,1)*temp(jm2,2)
     &                                    +temp(jp2,1)*temp(jp2,2))
     &                                +d1*(temp(jm1,1)*temp(jm1,2)
     &                                    +temp(jp1,1)*temp(jp1,2))
     &                                +d0*temp(l,1)*temp(l,2))
                f_rov =f_rov-filj*sigma_f*(d2*(temp(jm2,1)*temp(jm2,3)
     &                                    +temp(jp2,1)*temp(jp2,3))
     &                                +d1*(temp(jm1,1)*temp(jm1,3)
     &                                    +temp(jp1,1)*temp(jp1,3))
     &                                +d0*temp(l,1)*temp(l,3))
                f_row =f_row-filj*sigma_f*(d2*(temp(jm2,1)*temp(jm2,4)
     &                                    +temp(jp2,1)*temp(jp2,4))
     &                                +d1*(temp(jm1,1)*temp(jm1,4)
     &                                    +temp(jp1,1)*temp(jp1,4))
     &                                +d0*temp(l,1)*temp(l,4))

    !             ! selon k
    !             f_ro = f_ro-filk*sigma_f*(d2*(temp(km2,1)+temp(kp2,1))
    !  &                             +d1*(temp(km1,1)+temp(kp1,1))
    !  &                             +d0*temp(l,1))
    !             f_rou =f_rou-filk*sigma_f*(d2*(temp(km2,1)*temp(km2,2)
    !  &                                    +temp(kp2,1)*temp(kp2,2))
    !  &                                +d1*(temp(km1,1)*temp(km1,2)
    !  &                                    +temp(kp1,1)*temp(kp1,2))
    !  &                                +d0*temp(l,1)*temp(l,2))
    !             f_rov =f_rov-filk*sigma_f*(d2*(temp(km2,1)*temp(km2,3)
    !  &                                    +temp(kp2,1)*temp(kp2,3))
    !  &                                +d1*(temp(km1,1)*temp(km1,3)
    !  &                                    +temp(kp1,1)*temp(kp1,3))
    !  &                                +d0*temp(l,1)*temp(l,3))
    !             f_row =f_row-filk*sigma_f*(d2*(temp(km2,1)*temp(km2,4)
    !  &                                    +temp(kp2,1)*temp(kp2,4))
    !  &                                +d1*(temp(km1,1)*temp(km1,4)
    !  &                                    +temp(kp1,1)*temp(kp1,4))
    !  &                                +d0*temp(l,1)*temp(l,4))

                ! on met a jour les variables avec les grandeurs filtrees
                rop(l,1) = f_ro
                rop(l,2) = f_rou/f_ro
                rop(l,3) = f_rov/f_ro
                rop(l,4) = f_row/f_ro

              enddo
            enddo
           enddo


      end
