c***********************************************************************
c     $Date: 2010-11-04 13:25:50 +0100 (Thu, 04 Nov 2010) $
c     $Revision: 64 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine lesvist(ndo, param_int, param_real, neq_rot,depth,
     &                   ithread,nitrun,
     &                   ind_grad, ind_coe, ind_dm_zone,
     &                   xmut, rop, ti,tj,tk, vol, rot)

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
c_I    vent     : vitesses d entrainement aux facettes preced.
c
c     LOC
c_L    flu      : flux convectifs dans une direction de maillage
c
c     I/O
c_/    rot    : increment de la solution
c
c***********************************************************************
      implicit none

      INTEGER_E ndo, neq_rot, depth, ind_grad(6), ind_dm_zone(6),
     & ithread, nitrun, ind_coe(6), param_int(0:*)

      REAL_E rop(*), ti(*),tj(*),tk(*),vol(*), xmut(*), rot(*)

      REAL_E param_real(0:*)
C Var loc 
      INTEGER_E ind_flt(6), ind_extrap(6)
#include "FastS/param_solver.h"

          ind_flt(1) =max( ind_coe(1) , 0 )
          ind_flt(3) =max( ind_coe(3) , 0 )
          ind_flt(5) =max( ind_coe(5) , 0 )
          ind_flt(2) =min( ind_coe(2) , param_int(IJKV  ) +1 )
          ind_flt(4) =min( ind_coe(4) , param_int(IJKV+1) +1 )
          ind_flt(6) =min( ind_coe(6) , param_int(IJKV+2) +1 )

          ind_extrap(1) =max( ind_grad(1) , 1 )
          ind_extrap(3) =max( ind_grad(3) , 1 )
          ind_extrap(5) =max( ind_grad(5) , 1 )
          ind_extrap(2) =min( ind_grad(2) , param_int(IJKV  ) )
          ind_extrap(4) =min( ind_grad(4) , param_int(IJKV+1) )
          ind_extrap(6) =min( ind_grad(6) , param_int(IJKV+2) )

          !!IN: primitive var + metric
          !!out: rot(1)= rotx
          !!out: rot(2)= roty
          !!out: rot(3)= rotz
          !!out: xmut = (dudx**2 +dvdy**2 +dwdz**2 
          !!           + (dvdx+dudy)**2
          !!           + (dwdx+dudz)**2
          !!           + (dwdy+dvdz)**2)*vol*vol
          !!  
          !!  OUt valable sur une rangee de ghost
          call cp_rot_sij(ndo, param_int, neq_rot,
     &                    ind_flt,
     &                    xmut, rop, ti,tj,tk, vol, rot)
            
          !!  OUt : xmut = mulam + mu_sgs
          !!  OUt valable sur une rangee de ghost
          call mcmvi(ndo, param_int, param_real, neq_rot,ithread,nitrun,
     &                    ind_grad,
     &                    xmut, rop, rot)

          call extrap_mut(ndo, param_int, depth, 
     &                    ind_extrap, ind_dm_zone, xmut)
               
      end
