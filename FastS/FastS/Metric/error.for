c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 34 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine error(motc,noerr,itype)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Affichage centralise des erreurs
c
c     VAL
c
c     INP
c_I    motc : nom de la subroutine source d'erreur
c_I    noerr  : numero de l'erreur
c_I    itype  : type d'erreur (0:warning ; sinon:erreur et arret)
c
c***********************************************************************
      implicit none

#ifndef E_INTEGER_MPI_4
#  define      INTEGER_MPI integer*8
#else
#  define      INTEGER_MPI integer*4
#endif

      character *(*) motc
      INTEGER_E noerr,itype
c Var loc
      INTEGER_MPI ierr,ierr1
      INTEGER_E  iszmot,ithread,imp
c
c-----Numero de l'erreur et nom de la subroutine de provenance
c
      imp = 6
      iszmot=index(motc,"$",.true.)
      if(iszmot.eq.0) iszmot = 20

      write(imp,1000) noerr,motc(1:iszmot)
1000  format(' Message numero : ',i4,' - Subroutine : ',a10)
c
c
c-----Warning (0) ou Erreur (1)
c
      if (itype.eq.0) then
        write(imp,1010)
1010    format(' WARNING')
      else
        write(imp,1011)
1011    format(' ERROR')
      endif
c
c
c-----Messages relatifs a un sous dimensionnement de tableaux
c
      if (noerr.eq.1) then
        write(imp,9001)
      elseif (noerr.eq.2)  then
        write(imp,9002)
      elseif (noerr.eq.8) then
        write(imp,9008)
      elseif (noerr.eq.3) then
        write(imp,9003)
      elseif (noerr.eq.5) then
        write(imp,9005)
      elseif (noerr.eq.6) then
        write(imp,9006)
      elseif (noerr.eq.13) then
        write(imp,9013)
      elseif (noerr.eq.37) then
        write(imp,9037)
      elseif (noerr.eq.110) then
        write(imp,9110)
c
c
c-----Messages relatifs a une erreur de choix
c
      elseif (noerr.eq.4) then
        write(imp,9004)
      elseif (noerr.eq.7) then
        write(imp,9007)
      elseif (noerr.eq.16) then
        write(imp,9016)
      elseif (noerr.eq.18) then
        write(imp,9018)
      elseif (noerr.eq.19) then
        write(imp,9019)
      elseif (noerr.eq.20) then
        write(imp,9020)
      elseif (noerr.eq.21) then
        write(imp,9021)
      elseif (noerr.eq.22) then
        write(imp,9022)
      elseif (noerr.eq.23) then
        write(imp,9023)
      elseif (noerr.eq.25) then
        write(imp,9025)
      elseif (noerr.eq.26) then
        write(imp,9026)
      elseif (noerr.eq.27) then
        write(imp,9027)
      elseif (noerr.eq.28) then
        write(imp,9028)
      elseif (noerr.eq.30) then
        write(imp,9030)
      elseif (noerr.eq.31) then
        write(imp,9031)
      elseif (noerr.eq.33) then
        write(imp,9033)
      elseif (noerr.eq.34) then
        write(imp,9034)
      elseif (noerr.eq.35) then
        write(imp,9035)
      elseif (noerr.eq.38) then
        write(imp,9038)
      elseif (noerr.eq.39) then
        write(imp,9039)
c
c-----Messages relatifs a la coherence du schema numerique
c
      elseif (noerr.eq.9)  then
        write(imp,9009)
      elseif (noerr.eq.10) then
        write(imp,9010)
      elseif (noerr.eq.11) then
        write(imp,9011)
      elseif (noerr.eq.12) then
        write(imp,9012)
      elseif (noerr.eq.14) then
        write(imp,9014)
      elseif (noerr.eq.15) then
        write(imp,9015)
      elseif (noerr.eq.17) then
        write(imp,9017)
      elseif (noerr.eq.36) then
        write(imp,9036)
c
c-----Message relatif a une erreur de maintenace
c
      elseif (noerr.eq.70) then
        write(imp,9970)
      elseif (noerr.eq.71) then
        write(imp,9971)
      elseif (noerr.eq.90) then
        write(imp,9990)
c
c-----Erreur d'implementation de routine
c
      elseif (noerr.eq.61) then
        write(imp,9061)
c
c-----Erreur fatale demandant l'arret du programme
c
      elseif (noerr.eq.60) then
        write(imp,9060)
c
c-----Message relatif a un sous dimensionnement memoire
c
      elseif (noerr.eq.98) then
        write(imp,9998)
c
      elseif (noerr.eq.99) then
        write(imp,9998)
c
c
c-----Message relatif a une incapacite du code Magnus
c
      elseif (noerr.eq.200) then 
        write(imp,3000)
c
      elseif (noerr.eq.201) then
        write(imp,3001)
c
      elseif (noerr.eq.202) then
        write(imp,3002)
      
      elseif (noerr.eq.203) then
        write(imp,3003)

      elseif (noerr.eq.204) then
        write(imp,3004)

      elseif (noerr.eq.205) then
        write(imp,3005)

      elseif (noerr.eq.206) then
        write(imp,3006)
c
c  Erreurs relatives au schema differences finies
c
      elseif (noerr.eq.300) then
        write(imp,4300)
      elseif (noerr.eq.301) then
        write(imp,4301)
      elseif (noerr.eq.302) then
        write(imp,4302)
      elseif (noerr.eq.303) then
        write(imp,4303)
      elseif (noerr.eq.304) then
        write(imp,4304)
      elseif (noerr.eq.305) then
        write(imp,4305)
      elseif (noerr.eq.306) then
        write(imp,4306)
      elseif (noerr.eq.308) then
        write(imp,4308)

c
c  Erreurs relatives au terme source 
c
      elseif (noerr.eq.401) then
        write(imp,4401)
c
      elseif (noerr.eq.402) then
        write(imp,4402)

c
c  Erreurs relatives mpi
       elseif (noerr.eq.501) then
        write(imp,5001)
       elseif (noerr.eq.502) then
        write(imp,5002)
       elseif (noerr.eq.503) then
        write(imp,5003)
       elseif (noerr.eq.504) then
        write(imp,5004)

c  Erreurs relatives mpi
       elseif (noerr.eq.505) then
        write(imp,5005)
c
      else
        write(imp,'(a)') ' Pas de message pour cette erreur.'
      endif
c
c-----Warning (0) ou Erreur (1)
c
      if (itype.ne.0) then
           write(*,*)'ERROR: ithread=', motc(1:10)
           flush(6)
         stop
C         ithread  = 1
C!$       ithread  = OMP_get_thread_num() +1
C         if(error_count(ithread).eq.0) then
C           error_count(ithread) = 1
C           error_msg(ithread)   = motc(1:10)
C           write(*,*)'ERROR: ithread=',ithread, error_msg(ithread)
c         endif
      endif
c
c
c-----Formats d'ecriture
c
5001  format(' Numero de process trop grand.'/
     1       '  Redimensionner mxproc.'/)
5002  format(' Nbr de raccord mpi trop grand.'/
     1       '  Redimensionner mxracc.'/)
5003  format(' Nbr de zone laminaire trop grand.'/
     1       '  Redimensionner mx_zonelam dans SA.h.'/)
5004  format(' Maille au volume negatif...')
5005  format(' Nbr de zone IBC trop grand.'/
     1       '  Redimensionner mxibc dans com_srcterm.h.'/)
9001  format(' Numero de frontiere trop grand.'/
     1       '  Redimensionner mxbdr_no.'/)
9008  format(' Nombre de frontiere trop grand.'/
     1       '  Redimensionner mxbdr.'/)
9003  format(' Numero d''etat trop grand.'/
     1       '  Redimensionner  mxsta.'/)
9002  format(' Nombre de mot trop grand.'/
     1       '  Redimensionner  mxmot.'/)
9005  format(' Nombre d''options trop grand.'/
     1       '  Redimensionner  mxopt.'/)
9006  format(' Numero de domaine trop grand.'/
     1       '  Redimensionner  mxdom.'/)
9037  format(' Nombre de tableaux de la base de donnee trop grand.'/
     1       '  Redimensionner  mxdb.'/)
9013  format(' Nombre de lignes dans la boucle trop grand.'/
     1       '  Redimensionner  mxlin.'/)
9110  format(' Nombre de nacelles trop grand.'/
     1       '  Redimensionner  mxnac.'/)
c
9007  format(' Cette commande est inutile.'/'  Commande ignoree.'/)
c
9004  format(' Il faut, pour cette commande, selectionner un plan.'/
     1       ' Par defaut, le plan i2=i1 a ete choisi.'/)
c
9009  format(' Vous ne pouvez pas demander un calcul instationnaire'/
     1       ' et iterer plus d''une fois sur chaque plan.'/
     2       '  nitk a ete force a 1 .'/)
c
9010  format(' Vous ne pouvez pas faire un calcul instationnaire et'/
     1       ' utiliser l''ordre 0 .'/'  Mettre  ordre  a 1 ou 2 .'/)
9011  format(' Vous ne pouvez pas faire un calcul instationnaire et'/
     1       ' calculer le pas de temps localement.'/
     2       '  timestep a ete force a domain.'/)
9012  format(' Vous ne pouvez pas faire un calcul plan par plan sans'/
     1       ' avoir le meme nombre de plans dans chaque domaine.'/
     2       '  Modifier votre maillage.'/)
9014  format(' Il est inutile de faire plusieurs fois le meme calcul.'/
     1       '  nitp a ete force a 1'/)
9016  format(' Incoherence dans les entres/sorties pour les commandes'/
     1       ' read et EOF .'/
     2       ' Revoir les fichier de lecture de commandes.'/)
9018  format(' Le tableau que vous desirez liberer n''existe pas.'/
     1       '  La commande a ete ignoree.'/)
9019  format(' Incoherence dans les tailles de tableaux a liberer.'/)
9020  format(' La definition d''un etat de reference n''est pas ',
     1       ' implementee'/' pour un gaz reel.',/
     2       '  Cette commande a ete ignoree.'/)
9021  format(' Il faut selectionner un plan pour definir une ',
     1       ' condition limite.'/)
9022  format(' Cette condition limite n''est pas implementee pour',
     1       ' le type de fluide utilise.'/)
9023  format(' Il est conseille de preciser la condition limite en',/
     1       ' utilisant l''option (-i1, -i2, -j1, -j2, -k1 ou -k2).'/)
9025  format(' La taille d''un mot depasse 40 caracteres.'/)
9026  format(' Le nombre d''especes ne semble pas connue,'/
     1       ' il m''est impossible',
     2       ' de trouver la variable de l''energie.'/
     3       '  Option ignoree.'/)
9027  format(' La derniere commande a ete ignoree.'/)
9028  format(' Commande incompatible avec les precedentes.'/
     1       '  Verifier les dimensions des domaines.'/)
9030  format(' Type d''enregistrement inconnu.'/
     1       '  Verifier le fichier de lecture.'/)
9031  format(' Le fichier maillage ou flow n''a pas pu etre lu.'/
     1       '  Verifier le fichier de lecture.'/)
9033  format(' Cette commande necessite un domaine de redirection.'/
     1       '  Elle a ete ignoree.'/)
9034  format(' Vous n''avez pas defini le type d''iteration a',
     1       ' utiliser.'/
     2       '  L''iteration a ete forcee a  nitk .'/)
9035  format(' Cette declaration d''etat thermodynamique n''est pas'/
     1       ' implementee pour le type de gaz utilise.'/
     2       '  Modifier la declaration.'/)
9038  format(' L''option quick n''est pas implementee '/
     1       ' dans le cas de la modelisation k-epsilon.'/)
9039  format(' Seul le limiteur de van Albada est implemente '/
     1       ' dans le cas de la modelisation k-epsilon bi-espece.'/)
c
9060  format(' Cette erreur sera fatale.'/)
9061  format(' Cette routine n''est pas implementee.'/)
c
9015  format(' Il ne faut pas iterer plus d''une fois sur un plan'/
     1       ' numerique lorsque l''on utilise le schema a  2 pas ou'/
     2       ' si l''on est en multi-domaines.'/
     3       '  nitk a ete force a 1 et nitkp a la valeur de nitk .'/)
9036  format(' Il est plus efficace d''iterer sur nitk.'/
     1       '  Mettre nitk a ete force a nitkp et nitkp a 1 .'/)
9017  format(' L''option -update ne sert qu''en calcul domain.'/
     1       '  Option ignoree.'/)
c
9970  format(' Internal error : consult maintenance'/)
9971  format(' La condition d''entree par fichier n''est pas',
     1       ' completement implementee'/' dans cette version.'/
     2       '  Il est conseille d''utiliser a la place une condition',
     3       ' de raccord'/'  avec un domaine situe en amont.'/)
9990  format(' Fatal error : consult maintenance'/)
c
9998  format(' Fatal error : verify memory size')
c
3000  format('LU fm only'/)
3001  format('Unsteady: global and cstt dt >0 only'/)
3002  format('Steady: dt cstt>0 or global computation'/)
3003  format('constant dt only'/)
3004  format('Grid and metrics not updated: N instead of N+1')
3005  format('Grid and metrics updated: N+2 instead of N+1')
3006  format('Tableau pour le 2nd ordre en temps non cree')
4300  format('Schema valable uniquement a l''ordre 6')
4301  format('Rajouter des mailles fictives ou baisser l''ordre du'/ 
     1       ' schema.'/)
4302  format('Rajouter des mailles fictives ou baisser l''ordre du'/
     1       ' filtre.'/)
4303  format('Option 2d incompatible avec LES')
4304  format('Rajouter -2d apres df_centre (uniquement en DF), '/
     1       'ou mettre 2 plans minimum en k '/)
4305  format('Calcul de dro0 impossible avec l''option 2d/3d dans '/
     1       'le ''define fwsep'''/)
4306  format('Schema ECL avec ordre 4 uniquement')
4308  format('Option DF incompatible avec SA')

4401  format('Zone eponge : champ moyen non cree ; faire un create '/
     1       'statistiques + compute statistiques')
4402  format('Zone eponge : trop de zones, augmenter mxspg')

      end
