c***********************************************************************
c     $Date: 2010-06-14 09:57:46 +0200 (lun 14 jun 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine verrou(lok , ival, compteur)
c***********************************************************************
      implicit none

      INTEGER_E lok, ival
      integer*8 compteur

#include "parallelF.h"

C var loc
      INTEGER_E lok_tmp

      compteur = 0

      IF(ival.eq.1) then  !! verrou go

  100   continue

!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok
        if(lok_tmp.ne.0.and.lok_tmp.ne.2) then
           compteur = compteur + 1 
           goto 100
        endif
#ifdef _OPENMP4
!$OMP ATOMIC WRITE
#endif
        lok = 1
#ifdef _OPENMP4
!$OMP END ATOMIC 
#endif
!$OMP FLUSH (lok)

      ELSEIF(ival.eq.2) THEN    !! verrou wait scater

  200   continue
!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok
        if(lok_tmp.ne.1) then
          compteur = compteur + 1 
          goto 200
        endif
#ifdef _OPENMP4
!$OMP ATOMIC WRITE
#endif
        lok  = 2
#ifdef _OPENMP4
!$OMP END ATOMIC 
#endif
!$OMP FLUSH (lok)

      ELSE    !! verrou wait oldschool (ival.eq.3)


  300   continue
!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok
        if(lok_tmp.ne.1) then
          compteur = compteur + 1 
          goto 300
        endif

      ENDIF
      end

c***********************************************************************
c     $Date: 2010-06-14 09:57:46 +0200 (lun 14 jun 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine verrou_c(lok , ival)
c***********************************************************************
      implicit none

      INTEGER_E lok(1), ival

#include "parallelF.h"

C var loc
      INTEGER_E lok_tmp
      integer*8 compteur

      IF(ival.eq.1) then  !! verrou go

  100   continue

!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok(1)
        if(lok_tmp.ne.0.and.lok_tmp.ne.2) then
           goto 100
        endif
#ifdef _OPENMP4
!$OMP ATOMIC WRITE
#endif
        lok(1) = 1
#ifdef _OPENMP4
!$OMP END ATOMIC 
#endif
!$OMP FLUSH (lok)

      ELSEIF(ival.eq.2) THEN    !! verrou wait scater

  200   continue
!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok(1)
        if(lok_tmp.ne.1) then
          goto 200
        endif
c#ifdef _OPENMP4
c!$OMP ATOMIC WRITE
c#endif
c        lok(1)  = 2
c#ifdef _OPENMP4
c!$OMP END ATOMIC 
c#endif
c!$OMP FLUSH (lok)

      ELSEIF(ival.eq.3) THEN   !! verrou wait oldschool (ival.eq.3)


  300   continue
!$OMP FLUSH (lok)
#ifdef _OPENMP4
!$OMP ATOMIC READ
#endif
        lok_tmp = lok(1)
        if(lok_tmp.ne.1) then
          goto 300
        endif

      ELSE  !! init val

#ifdef _OPENMP4
!$OMP ATOMIC WRITE
#endif
        lok(1) = 0
#ifdef _OPENMP4
!$OMP END ATOMIC 
#endif
!$OMP FLUSH (lok)

      ENDIF
      end

