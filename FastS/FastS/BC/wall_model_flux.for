c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine wall_model_flux(ndom, ithread, param_int, param_real, 
     &                  ind_dm_zone, ind_sdm, nitcfg, nitrun, cycl,
     &                  psi, wig, stat_wig, rop, drodm, x,y,z,
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
     & venti(*),ventj(*),ventk(*), wig(*),stat_wig(*),x(*),y(*),z(*),
     & ti_df(*),tj_df(*),tk_df(*),vol_df(*), psi(*), param_real(0:*)

C     Var loc
      INTEGER_E ndf,lf,npass,incijk, nb_bc, iptdata, flag_correct_flu
      INTEGER_E bc_type, pt_bc,pt_bcs, idir,nbdata,lskip,size_data
      INTEGER_E ind_CL(6), ind_CL119(6), nitcfg, cycl
      INTEGER_E sampling, nitrun, iflow,lcompute
      REAL_E mobile_coef
 
      ! pas de debordement coin pour calcul des flu
      npass =0
      pt_bcs = param_int(PT_BC)
      nb_bc = param_int(pt_bcs)
      do 100 ndf = 0, nb_bc-1
         
        pt_bc  = param_int(pt_bcs + 1 + ndf)

        bc_type= param_int(pt_bc + BC_TYPE)

        !write(*,*)'ndf',bc_type,ndf,ndom, nb_bc-1

        !on corrige wall et symetrie pour le moment
        !if(bc_type.ne.BCWALLMODEL.and.bc_type.ne.BCWALLEXCHANGE) goto 100
        if(bc_type.ne.BCWALLMODEL) goto 100
        !if(bc_type.ne.BCWALLEXCHANGE) goto 100
        !write(*,*)'ndf',bc_type,ndf,ndom, nb_bc-1

        idir   = param_int( pt_bc + BC_IDIR)
        nbdata = param_int( pt_bc + BC_NBDATA)
        size_data = param_int(ndom) + pt_bc + BC_NBDATA + 1
        iptdata = 0
        if (nbdata.ne.0) iptdata = param_int( pt_bcs + 1 + ndf + nb_bc)

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

         cycl=param_int(NSSITER)/param_int(LEVEL)
c        write(*,*)'avt correction',ndom,nitcfg
c     & , param_int(LEVEL)

        ! ---- Correction of the boundary face fluxes

        if(param_int(EXPLOC).ne.0.and.mod(nitcfg,cycl).ne.1) goto 100 ! pas de calcul a cett ssiter en dtlocal

        ! on retranche flux convectif ET visqueux
        iflow = 2
        if(param_int(KFLUDOM).eq.2) then

          call corr_flusenseur_select(ndom, ithread, idir, iflow,
     &                        param_int, param_real,
     &                        ind_CL,
     &                        rop, drodm  , wig,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol, xmut)

        elseif(param_int(KFLUDOM).eq.1) then

          call corr_fluausm_select(ndom, ithread, idir, iflow,
     &                        param_int, param_real,
     &                        ind_CL,
     &                        rop, drodm  , wig,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol, xmut)

        elseif(param_int(KFLUDOM).eq.5) then

          call corr_fluroe_select(ndom, ithread, idir, iflow,
     &                        param_int, param_real,
     &                        ind_CL,
     &                        rop, drodm  , wig,
     &                        venti, ventj, ventk,
     &                        ti,tj,tk,vol, xmut)

        else
           call error('wallmodel_flu$',70,1)
        endif


#include "FastS/BC/CL_wallmodel_flu.for"

 100  continue

      end
