c***********************************************************************
c     $Date: 2010-06-14 09:57:46 +0200 (lun 14 jun 2010) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine flush_integer(isize,tab)
c***********************************************************************
      implicit none

      INTEGER_E isize, tab(isize)
!$OMP FLUSH (tab)
      end
