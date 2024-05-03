c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine bfl3(ndom, ithread, param_int, param_real, 
     &                  ind_dm_zone, ind_sdm,
     &                  psi, wig, stat_wig, rop, drodm,
     &                  ti,ti_df,tj,tj_df,tk,tk_df, vol,vol_df,
     &                  venti, ventj, ventk, xmut)
c***********************************************************************
c     _U   USER : PECHIER
c     
c     ACT
c     _A    Appel des routines de calcul des flux frontieres aux interfaces
c     _A    du domaine.
c     
c     VAL
c     _V    processeur domaine
c     _V    unsteady/steady
c     
c     OUT
c     _O    drodm
C     
c     CONDITION LIMITES PAR MODIFICATION DE FLUX FRONTIERE
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ithread,ind_dm_zone(6),ind_sdm(6), param_int(0:*)

      REAL_E rop(*),drodm(*),xmut(*), ti(*),tj(*),tk(*),vol(*),
     & venti(*),ventj(*),ventk(*), wig(*),stat_wig(*),
     & ti_df(*),tj_df(*),tk_df(*),vol_df(*), psi(*), param_real(0:*)

C     Var loc
      INTEGER_E ndf,lf,npass,incijk, nb_bc, iptdata, flag_correct_flu
      INTEGER_E bc_type, pt_bc,pt_bcs, idir,nbdata,lskip, iflow
      INTEGER_E ind_CL(6), ind_CL119(6)
      REAL_E mobile_coef
 
      ! pas de debordement coin pour calcul des flu
      npass =0
      pt_bcs = param_int(PT_BC)
      nb_bc = param_int(pt_bcs)
      !write(*,*)'bc', nb_bc
      do 100 ndf = 0, nb_bc-1
         
        pt_bc  = param_int(pt_bcs + 1 + ndf)

        bc_type= param_int(pt_bc + BC_TYPE)

        !write(*,*)'ndf',bc_type,ndf

        !on corrige wall et symetrie pour le moment
         If(    bc_type.ne.BCWALLINVISCID.and.bc_type.ne.BCSYMMETRYPLANE
     &     .and.bc_type.ne.BCWALL        .and.bc_type.ne.BCWALLVISCOUS
     &     .and.bc_type.ne.BCWALLVISCOUS_ISOT_FICH
     &     .and.bc_type.ne.BCWALLVISCOUS_TRANSITION)
C     &     .and.bc_type.ne.BCINFLOW.and.bc_type.ne.BCINJ1) 
     &   goto 100

        idir   = param_int( pt_bc + BC_IDIR)
        nbdata = param_int( pt_bc + BC_NBDATA)


       !calcul intersection sous domaine et fenetre CL
        call indice_cl_sdm(idir,npass,lskip,param_int( pt_bc+BC_TYPE),
     &                     param_int(NIJK+3), param_int(NIJK+4),
     &                     param_int(IJKV)  , param_int(pt_bc + BC_FEN),
     &                     ind_dm_zone, ind_sdm,                !IN
     &                     ind_CL, ind_CL119);                  !OUT



         if(lskip.ne.0) goto 100 !la fenetre sous-domaine n'intersecte pas la Cond Limite no(lf)

         if(idir.eq.1) then
           ind_CL(1)=1
           ind_CL(2)=1
         elseif(idir.eq.2) then
           ind_CL(2)=ind_CL(1)
         elseif(idir.eq.3) then
           ind_CL(3)=1
           ind_CL(4)=1
         elseif(idir.eq.4) then
           ind_CL(4)=ind_CL(3)
         elseif(idir.eq.5) then
           ind_CL(5)=1
           ind_CL(6)=1
         else
           ind_CL(6)=ind_CL(5)
         endif

         if (param_int(KFLUDOM).eq.5) then

           !!  attention flux convectif 6eme variable SA pas retranché
           !! donc pas de garantie de flux(6) =0 a la paroi en Roe
           iflow =1 
           
           call corr_fluroe_select(ndom, ithread, idir, iflow,
     &                        param_int, param_real,
     &                        ind_CL, 
     &                        rop, drodm  , wig,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol, xmut)
         else
           continue
         endif


#include "FastS/BC/CL_correction_flu.for"

        !c.....Traitement de non reflexion 
c         If (lnrfront(lf)) Then 
c 
c          call bflnr(ndom,idir,lf, 
c     &               neq,neq_mtr,neq_vent,ndimdx, 
c     &               inc2,inc3,
c     &               ind_CL,
c     &               rop, 
c     &               flu, 
c     &               tijk_df,    !metrique VF ou DF 
c     &               ventijk)
c          Endif

 100  continue
c

      end
