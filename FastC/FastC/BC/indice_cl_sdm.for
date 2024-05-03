c***********************************************************************
c     $Date: 2011-07-28 18:43:32 +0200 (jeu 28 jui 2011) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine indice_cl_sdm(idir, npass, lskip, type_bc,
     &                       ific_ij_loc, ific_k_loc,
     &                       ijkv, ind_fen, ind_dm, ind_sdm,
     &                       ind_CL, ind_CL119)
c***********************************************************************
c_U  USER : PECHIER
c
c     ACTc  - Paroi adiabatique                  : bgad3
c_A    calcul des dimension des fenetre de condition limite sur le sous domaine courant
c
c     VAL
c_V    Valable pour le processeur domaine
c
c     INP
c_I    lpass: vrai si on traite les raccord  (2eme appel a init_var_dom dans int_var_omp)
c_I    npass: nombre de passage totale dans les CL
c_I          : 2 passage necessaire pour remplir correctement les coins des raccord match
c_I    ipass: indice de boucle sur les passage dans les CL. Ipass varie de 1 a npass pour les extract et vaut 1 pour
c_I           les comput

c
c     OUT
c_O   ind_CL, lskip
c=======================================================================
      implicit none


      INTEGER_E lskip, idir,npass,ific_ij_loc,ific_k_loc,
     &  ijkv(3),ind_fen(6), ind_sdm(6),ind_CL(6),ind_CL119(6), ind_dm(6)
      INTEGER_E type_bc

C Var loc
      INTEGER_E ic1


c      write(*,*)'ijkv', ijkv
c      write(*,*)'fen',ind_fen
c      write(*,*)'sdm',ind_sdm
c      write(*,*)'dm',ind_dm

      ind_CL(1) = 1
      ind_CL(2) = 0
      ind_CL(3) = 1
      ind_CL(4) = 0
      ind_CL(5) = 1
      ind_CL(6) = 0

      IF(     (ind_sdm(1).le.ind_fen(2)).and.(ind_sdm(2).ge.ind_fen(1))
     &   .and.(ind_sdm(3).le.ind_fen(4)).and.(ind_sdm(4).ge.ind_fen(3)) 
     &   .and.(ind_sdm(5).le.ind_fen(6)).and.(ind_sdm(6).ge.ind_fen(5)))
     &    THEN

       ind_CL(1)=max(ind_fen(1),ind_sdm(1))
       ind_CL(2)=min(ind_fen(2),ind_sdm(2))
       ind_CL(3)=max(ind_fen(3),ind_sdm(3))
       ind_CL(4)=min(ind_fen(4),ind_sdm(4))
       ind_CL(5)=max(ind_fen(5),ind_sdm(5))
       ind_CL(6)=min(ind_fen(6),ind_sdm(6))

        if(idir.eq.1) then
           ind_CL(2)    = ind_CL(1)- 1
           ind_CL(1)    = ind_CL(1)- ific_ij_loc
        elseif(idir.eq.2) then
           ind_CL(2)    = ind_CL(1)+ ific_ij_loc
           ind_CL(1)    = ind_CL(1)+ 1
        elseif(idir.eq.3) then
           ind_CL(4)    = ind_CL(3)- 1
           ind_CL(3)    = ind_CL(3)- ific_ij_loc
        elseif(idir.eq.4) then
             ind_CL(4)    = ind_CL(3)+ ific_ij_loc
             ind_CL(3)    = ind_CL(3)+ 1
        elseif(idir.eq.5) then
             ind_CL(6)    = ind_CL(5)- 1
             ind_CL(5)    = ind_CL(5)- ific_k_loc
        else 
             ind_CL(6)    = ind_CL(5)+ ific_k_loc
             ind_CL(5)    = ind_CL(5)+ 1
        endif


       ind_CL119(1) =  ind_CL(1)
       ind_CL119(2) =  ind_CL(2)
       ind_CL119(3) =  ind_CL(3)
       ind_CL119(4) =  ind_CL(4)
       ind_CL119(5) =  ind_CL(5)
       ind_CL119(6) =  ind_CL(6)

       lskip = 0

       !on deborde dans les coin ssi necessaire:
       !  - 2 passes et bord domaine
       !  - bord sous domaine (LU_loc) si npass=1 ou 2

       if(npass.eq.0) goto 1000 !zero debordement 

        ic1 = 1
C Nearmatch
      if(type_bc.eq.9.and.npass.eq.2) ic1= 0
C BCWallViscous_isot_fich
      if(type_bc.eq.7.and.npass.eq.2)ic1=0
C BC_fich
      if(type_bc.eq.8.and.npass.eq.2) ic1=0
C BC_farfield (cart normale maille coin = bofbof)
      !if(type_bc.eq.1.and.npass.eq.2) ic1=0

       !   if(icle.eq.-734.and.npass.eq.2) ic1= 0
       !   if(icle.eq. 230.and.npass.eq.2) ic1= 0  ! test sortie routine pour empecher depassement
       !   if(icle.eq.-399.and.npass.eq.2) ic1= 0  ! test sortie routine pour empecher depassement

       IF(idir.eq.1.or.idir.eq.2) then 

           if(  ( ind_CL(3).eq.ind_dm(3) ) .and.
     &          (    (npass.eq.1.and.ind_CL(3).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(3).eq.1) ) ) then

               ind_CL(3)   = ind_CL(3)   -ific_ij_loc*ic1
               ind_CL119(3)= ind_CL119(3)-ific_ij_loc
           endif
           if(  ( ind_CL(4).eq.ind_dm(4) ) .and.
     &          (    (npass.eq.1.and.ind_CL(4).ne.ijkv(2))
     &           .or.(npass.eq.2.and.ind_CL(4).eq.ijkv(2)) ) ) then

               ind_CL(4)   = ind_CL(4)   +ific_ij_loc*ic1
               ind_CL119(4)= ind_CL119(4)+ific_ij_loc
           endif
           if(  ( ind_CL(5).eq.ind_dm(5) ) .and.
     &          (    (npass.eq.1.and.ind_CL(5).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(5).eq.1) ) ) then

               ind_CL(5)   = ind_CL(5)   -ific_k_loc*ic1
               ind_CL119(5)= ind_CL119(5)-ific_k_loc
           endif
           if(  ( ind_CL(6).eq.ind_dm(6) ) .and.
     &          (    (npass.eq.1.and.ind_CL(6).ne.ijkv(3))
     &           .or.(npass.eq.2.and.ind_CL(6).eq.ijkv(3)) ) ) then

               ind_CL(6)   = ind_CL(6)   +ific_k_loc*ic1
               ind_CL119(6)= ind_CL119(6)+ific_k_loc
           endif

          !empeche debordement si sous-domaine et frontiere fichier multi-fenetre
C BC_fich ou BCWallViscous_isot_fich
          if(type_bc.eq.8.or.
     &       type_bc.eq.7) then
           ind_CL(3) = max(ind_CL(3),ind_fen(3))
           ind_CL(4) = min(ind_CL(4),ind_fen(4))
           ind_CL(5) = max(ind_CL(5),ind_fen(5))
           ind_CL(6) = min(ind_CL(6),ind_fen(6))
          endif

          !if(idir.eq.1) then
          !   ind_CL(2)    = ind_CL(1)- 1
          !   ind_CL(1)    = ind_CL(1)- ific_ij_loc
          !   ind_CL119(2) = ind_CL(2)
          !   ind_CL119(1) = ind_CL(1)
          !else 
          !   ind_CL(2)    = ind_CL(1)+ ific_ij_loc
          !   ind_CL(1)    = ind_CL(1)+ 1
          !   ind_CL119(2) = ind_CL(2)
          !   ind_CL119(1) = ind_CL(1)
          !endif

       ELSEIF(idir.eq.3.or.idir.eq.4) then 

           if(  ( ind_CL(1).eq.ind_dm(1) ) .and.
     &          (    (npass.eq.1.and.ind_CL(1).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(1).eq.1) ) ) then

               ind_CL(1)   = ind_CL(1)   -ific_ij_loc*ic1
               ind_CL119(1)= ind_CL119(1)-ific_ij_loc
           endif
           if(  ( ind_CL(2).eq.ind_dm(2) ) .and.
     &          (    (npass.eq.1.and.ind_CL(2).ne.ijkv(1))
     &           .or.(npass.eq.2.and.ind_CL(2).eq.ijkv(1)) ) ) then

               ind_CL(2)   = ind_CL(2)   +ific_ij_loc*ic1
               ind_CL119(2)= ind_CL119(2)+ific_ij_loc
           endif
           if(  ( ind_CL(5).eq.ind_dm(5) ) .and.
     &          (    (npass.eq.1.and.ind_CL(5).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(5).eq.1) ) ) then

               ind_CL(5)   = ind_CL(5)   -ific_k_loc*ic1
               ind_CL119(5)= ind_CL119(5)-ific_k_loc
           endif
           if(  ( ind_CL(6).eq.ind_dm(6) ) .and.
     &          (    (npass.eq.1.and.ind_CL(6).ne.ijkv(3))
     &           .or.(npass.eq.2.and.ind_CL(6).eq.ijkv(3)) ) ) then

               ind_CL(6)   = ind_CL(6)   +ific_k_loc*ic1
               ind_CL119(6)= ind_CL119(6)+ific_k_loc
           endif

          !empeche debordement si sous-domaine et frontiere fichier multi-fenetre
C BC_fich ou BCWallViscous_isot_fich
          if(type_bc.eq.8.or.
     &       type_bc.eq.7) then
           ind_CL(1) = max(ind_CL(1),ind_fen(1))
           ind_CL(2) = min(ind_CL(2),ind_fen(2))
           ind_CL(5) = max(ind_CL(5),ind_fen(5))
           ind_CL(6) = min(ind_CL(6),ind_fen(6))
          endif

          !if(idir.eq.3) then
          !   ind_CL(4)    = ind_CL(3)- 1
          !   ind_CL(3)    = ind_CL(3)- ific_ij_loc
          !   ind_CL119(4) = ind_CL(4)
          !   ind_CL119(3) = ind_CL(3)
          !else 
          !   ind_CL(4)    = ind_CL(3)+ ific_ij_loc
          !   ind_CL(3)    = ind_CL(3)+ 1
          !   ind_CL119(4) = ind_CL(4)
          !   ind_CL119(3) = ind_CL(3)
          !endif

       ELSE
           if(  ( ind_CL(1).eq.ind_dm(1) ) .and.
     &          (    (npass.eq.1.and.ind_CL(1).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(1).eq.1) ) ) then

               ind_CL(1)   = ind_CL(1)   -ific_ij_loc*ic1
               ind_CL119(1)= ind_CL119(1)-ific_ij_loc
           endif
           if(  ( ind_CL(2).eq.ind_dm(2) ) .and.
     &          (    (npass.eq.1.and.ind_CL(2).ne.ijkv(1))
     &           .or.(npass.eq.2.and.ind_CL(2).eq.ijkv(1)) ) ) then

               ind_CL(2)   = ind_CL(2)   +ific_ij_loc*ic1
               ind_CL119(2)= ind_CL119(2)+ific_ij_loc
           endif
           if(  ( ind_CL(3).eq.ind_dm(3) ) .and.
     &          (    (npass.eq.1.and.ind_CL(3).ne.1)
     &           .or.(npass.eq.2.and.ind_CL(3).eq.1) ) ) then

               ind_CL(3)   = ind_CL(3)   -ific_ij_loc*ic1
               ind_CL119(3)= ind_CL119(3)-ific_ij_loc
           endif
           if(  ( ind_CL(4).eq.ind_dm(4) ) .and.
     &          (    (npass.eq.1.and.ind_CL(4).ne.ijkv(2))
     &           .or.(npass.eq.2.and.ind_CL(4).eq.ijkv(2)) ) ) then

               ind_CL(4)   = ind_CL(4)   +ific_ij_loc*ic1
               ind_CL119(4)= ind_CL119(4)+ific_ij_loc
           endif

          !empeche debordement si sous-domaine et frontiere fichier multi-fenetre
C BC_fich ou BCWallViscous_isot_fich
          if(type_bc.eq.7.or.
     &       type_bc.eq.8) then
           ind_CL(1) = max(ind_CL(1),ind_fen(1))
           ind_CL(2) = min(ind_CL(2),ind_fen(2))
           ind_CL(3) = max(ind_CL(3),ind_fen(3))
           ind_CL(4) = min(ind_CL(4),ind_fen(4))
          endif

          !if(idir.eq.5) then
          !   ind_CL(6)    = ind_CL(5)- 1
          !   ind_CL(5)    = ind_CL(5)- ific_k_loc
          !   ind_CL119(6) = ind_CL(6)
          !   ind_CL119(5) = ind_CL(5)
          !else 
          !   ind_CL(6)    = ind_CL(5)+ ific_k_loc
          !   ind_CL(5)    = ind_CL(5)+ 1
          !   ind_CL119(6) = ind_CL(6)
          !   ind_CL119(5) = ind_CL(5)
          !endif

       endif!idir


1000  continue

      ELSE
          lskip = 1 !la fenetre sous-domaine n'intersecte pas la Cond Limite 
      ENDIF
 
      end
