c***********************************************************************
c     $Date: 2014-03-19 20:08:08 +0100 (mer. 19 mars 2014) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine correct_flux(ndom, ithread, param_int, param_real, 
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
      INTEGER_E ind_CL(6), ind_CL119(6), inc_bc(3), nitcfg, cycl
      INTEGER_E sampling, nitrun, iflow,lapply, iptfen
      REAL_E mobile_coef
 
      cycl=param_int(NSSITER)/param_int(LEVEL)
      if(param_int(EXPLOC).ne.0.and.mod(nitcfg,cycl).ne.1) return ! pas de calcul a cette ssiter en dtlocal


      ! pas de debordement coin pour calcul des flu
      npass =0
      pt_bcs = param_int(PT_BC)
      nb_bc = param_int(pt_bcs)
      do 100 ndf = 0, nb_bc-1
         
        pt_bc  = param_int(pt_bcs + 1 + ndf)

        bc_type= param_int(pt_bc + BC_TYPE)

        !On verifie si la CL necessite un flux particulier
        lapply = 0
        if(bc_type.eq.BCWALLMODEL                   ) then
          lapply=1
        elseif(bc_type.eq.BCFluxOctreeF ) then
          lapply=2
        elseif(bc_type.eq.BCFluxOctreeC.and.nitcfg.gt.1) then
          lapply=3
        endif

        If( (   bc_type.eq.BCWALLINVISCID.or.bc_type.eq.BCSYMMETRYPLANE
     &      .or.bc_type.eq.BCWALL        .or.bc_type.eq.BCWALLVISCOUS
     &      .or.bc_type.eq.BCWALLVISCOUS_ISOT_FICH
     &      .or.bc_type.eq.BCWALLVISCOUS_TRANSITION
     &      ) .and.param_int(KFLUDOM).eq.5)      lapply = 4

        if(lapply.eq.0) goto 100 !correction de flux pas necessaire

        idir   = param_int( pt_bc + BC_IDIR)
        nbdata = param_int( pt_bc + BC_NBDATA)
        size_data = param_int( pt_bc + BC_NBDATA + 1)/param_int(NEQ)
        !!size_data = param_int( pt_bc + BC_NBDATA + !1)/(param_int(NEQ)*2)  ! tentative pour envoie pente
        iptdata = 0
        if (nbdata.ne.0) iptdata = param_int( pt_bcs + 1 + ndf + nb_bc)

        !calcul intersection sous domaine et fenetre CL
        call indice_cl_sdm(idir,npass,lskip,param_int( pt_bc+BC_TYPE),
     &                     param_int(NIJK+3), param_int(NIJK+4),
     &                     param_int(IJKV)  , param_int(pt_bc + BC_FEN),
     &                     ind_dm_zone, ind_sdm,                !IN
     &                     ind_CL, ind_CL119);                  !OUT

c      if(ndom.eq.48) write(*,'(a,i2,a,i3,a,i3)')'skip',lskip,
c     & ' ndf=',ndf,' idir=', idir
c      if(ndom.eq.48) write(*,'(a,6i4)')'sdm',ind_sdm(1:6)
c      if(ndom.eq.48) write(*,'(a,i4,a,i3)')'ipdata', iptdata,
c     & 'size', size_data
c      if(ndom.eq.48) write(*,'(a,6i4)')'fen',
c     & param_int(pt_bc + BC_FEN: pt_bc + BC_FEN+5)

        if(lskip.ne.0) goto 100 !la fenetre sous-domaine n'intersecte pas la Cond Limite no(lf)

        iptfen = pt_bc + BC_FEN
        if(idir.le.2) then
           inc_bc(1)= param_int(iptfen +3) - param_int(iptfen +2) +1 !nombre element de la fenetre dans la direction J
           inc_bc(2)= param_int(iptfen +2) ! debut indice j
           inc_bc(3)= param_int(iptfen +4) ! debut indice k
        elseif(idir.le.4) then
           inc_bc(1)= param_int(iptfen +1) - param_int(iptfen ) +1 !nombre element de la fenetre dans la direction I
           inc_bc(2)= param_int(iptfen   ) ! debut indice i
           inc_bc(3)= param_int(iptfen +4) ! debut indice k
        else
           inc_bc(1)= param_int(iptfen +1) - param_int(iptfen ) +1 !nombre element de la fenetre dans la direction I
           inc_bc(2)= param_int(iptfen   ) ! debut indice i
           inc_bc(3)= param_int(iptfen +2) ! debut indice j
        endif
        
        !modifie range normal sauf pour flux conservatif
        If(lapply.ne.2.and.lapply.ne.3) Then
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
        Endif

        ! on retranche ou stocke flux convectif ET visqueux
        iflow = 2               ! on corrige juste les flux laminaires
        if(lapply.eq.4) iflow=1 ! on corrige juste les flux euler pour paroi +Roe 


c      write(*,'(a,i4,a,4i4,a,i3,a,i3,a,i3,a,i3)')
c     & 'nd=',ndom,' indCL=',ind_CL(1:4),' lapply=',
c     & lapply,' nstep=',nitcfg,' ndf=',ndf, ' idir=', idir
      !if(ndom.eq.0) write(*,*)'indbc', inc_bc(1:3)

        if(param_int(KFLUDOM).eq.2) then
 
          call corr_flusenseur_select(ndom, ithread, idir, bc_type,
     &                          size_data, iflow,
     &                          param_int, param_real,
     &                          ind_CL, inc_bc,
     &                          rop, drodm  , wig, param_real(iptdata),
     &                          venti, ventj, ventk,
     &                          ti,tj,tk,vol, xmut)

        elseif(param_int(KFLUDOM).eq.1) then

          call corr_fluausm_select(ndom, ithread, idir, bc_type,
     &                          size_data, iflow,
     &                          param_int, param_real,
     &                          ind_CL, inc_bc,
     &                          rop, drodm  , wig, param_real(iptdata),
     &                          venti, ventj, ventk,
     &                          ti,tj,tk,vol, xmut)

        elseif(param_int(KFLUDOM).eq.5) then


          call corr_fluroe_select(ndom, ithread, idir, bc_type,
     &                          size_data, iflow,
     &                          param_int, param_real,
     &                          ind_CL, inc_bc,
     &                          rop, drodm  , wig, param_real(iptdata),
     &                          venti, ventj, ventk,
     &                          ti,tj,tk,vol, xmut)

        else
             call error('correct_flu$',70,1)
        endif

     
      !stokage et correction flux fait dans dans fct corr... pour flux
      !conservatif; Reste wall(roe) et wallmodel a faire
      if(lapply.eq.1) then   
#include "FastS/BC/CL_wallmodel_flu.for"
      elseif(lapply.eq.4) then
#include "FastS/BC/CL_correction_flu.for"
      endif
      

 100  continue

      end
