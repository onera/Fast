c***********************************************************************
c     $Date: 2016-03-01 20:08:08 +0100 (mar. 1 mars 2016) $
c     $Revision : ?? $
c     $Author: DidierBlaise $
c***********************************************************************
      subroutine conv2pytree(cvg_ptr, nitrun, neq, LastRec, zone_name,
     &   size_name, lft, nrec, nd, Itnum, Res_L2, Res_oo)
c
c     ACT
c_A     Stockage des residus moyens du processus sous iteratif (L_2 et oo)
c
c     INP
c_I    cvg_ptr, nitrun, neq, LastRec, zone_name, size_name, lft
c
c     OUT
c_O   Itnum, Res_L2, Res_oo
c***********************************************************************
c UTILISATION
c                    display ss_iteration +ndoms
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      implicit none
      INTEGER_E nrec, nd
      INTEGER_E nitrun, neq, Itnum(nrec), LastRec, size_name, lft
      REAL_E cvg_ptr(2*neq), Res_L2(neq*nrec), Res_oo(neq*nrec)
      Character(len=size_name)zone_name 
      Character(len=40)zname
      Character(len=132)vart
c Var loc
      INTEGER_E i,ne, Irec, ifile,ilen,j
      INTEGER_E isize
C     STRUCTURE /recconv/
        INTEGER_E Itn(nrec)
        REAL_E ResL2(8,nrec)
        REAL_E Resoo(8,nrec)
C     END STRUCTURE
C     RECORD /recconv/ page(nrec)
C     
c
c     ACT
c_A     Stockage des residus moyens du processus sous iteratif (L_2)
c
c     INP
c_I    cvg_ptr, nit
C
      IF(LastRec .EQ . (nrec + 1) .OR. lft.LT.0)THEN
        DO I=1,LastRec
        Irec = (I-1)*neq
        Itn(I) = Itnum(I)
          DO J=1,neq
          ResL2(J,I) = Res_L2(Irec+J)
          Resoo(J,I) = Res_oo(Irec+J)
          ENDDO
        ENDDO
C
        WRITE(zname,'(40A)')(char(32),I=1,40)
        ilen= min(size_name,40)
        zname(1:ilen)=zone_name(1:ilen)
C
      IFILE=1955
      IF(lft.EQ.-1)THEN
       OPEN(UNIT=IFILE,FILE="residuals.bin",
     &      FORM="unformatted",POSITION="append")
        WRITE(IFILE)zname, neq, LastRec, 
     &  ( Itn(I),(ResL2(J,I),Resoo(J,I),J=1,neq),I=1,LastRec)
        CLOSE(IFILE)
        LastRec=0
      ELSE
        OPEN(UNIT=IFILE,FILE="residuals.dat",POSITION="append")
      vart='VARIABLES="It" "Ro_L2" "RoU_L2" "RoV_L2" "RoW_L2" "RoE_L2"'
          IF(neq.eq.6)THEN
          vart=TRIM(vart)//' "Nut_L2"'
          ENDIF
      vart=TRIM(vart)//' "Ro_oo" "RoU_oo" "RoV_oo" "RoW_oo" "RoE_oo"'
          IF(neq.eq.6)THEN
          vart=TRIM(vart)//' "Nut_oo"'
          ENDIF
          if(nd.eq.0)WRITE(IFILE,'(A)')TRIM(vart)
      WRITE(IFILE,'(A)')'ZONE T="'//zname(1:ilen)//'"'
        DO I=1,LastRec
        WRITE(IFILE,'(I7,16(1X,E15.8))')
     &   Itn(I),(ResL2(J,I),Resoo(J,I),J=1,neq)
        ENDDO
      CLOSE(IFILE)
      ENDIF
      ENDIF

      IF( lft.GE.0)THEN
      Irec = (neq-1)*LastRec
      LastRec = LastRec + 1
      Itnum(LastRec) = nitrun
        DO ne=1,neq
        Irec = (LastRec-1)*neq + ne
        Res_L2( Irec ) = cvg_ptr(ne)
        Res_oo( Irec ) = cvg_ptr(neq + ne)
        ENDDO
      ENDIF
C
C
C        F90 Version
C       DO I=1,LastRec
C       Irec = (I-1)*neq
C       page(I).Itn = Itnum(I)
C       page(I).ResL2(1:neq) = Res_L2(Irec+1:Irec+neq)
C       page(I).Resoo(1:neq) = Res_oo(Irec+1:Irec+neq)
C       ENDDO
C
C       write(zname,'(40A)')(char(32),I=1,40)
C       ilen= min(size_name,40)
C       zname(1:ilen)=zone_name(1:ilen)
C       print *, zname
C
C       IFILE=1955
C       OPEN(UNIT=IFILE,FILE="residuals.bin",
C    &       FORM="unformatted",POSITION="append")     
C       WRITE(IFILE)zname, neq, LastRec, page(1:LastRec)
C       CLOSE(IFILE)
C       LastRec=0
C     ENDIF
C
      end
