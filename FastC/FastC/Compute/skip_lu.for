c***********************************************************************
c     $Date: 2010-06-30 10:59:35 +0200 (Wed, 30 Jun 2010) $
c     $Revision: 34 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine skip_lu(nitmax, mx_ssiter, iskip_lu)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  igetadr.sf/2.0/e1/i1
c
c_DC  DATE_C : Sep  8 1993 -- AUTEUR : 00000000
c
c     HISTORIQUE
c_H    Iteration     1 ___ 1994-07-20 16:38:03  mletou
c_H    Creation

c
c     ACT
c_A    Cette fonction renvoie l'adresse du tableau numero 'no'.
c
c     VAL
c
c     INP
c_I    no      : numero du tableau dans la Base de Donnee
c
c     OUT
c_O    igetadr = adresse du tableau numero 'no'.
c***********************************************************************
      implicit none

      INTEGER_E nitmax, mx_ssiter,  iskip_lu(mx_ssiter, 2)

C Var loc
      INTEGER_E nq

      nq = 2   
      if(  nitmax .eq.  mx_ssiter) nq = 1   
 
      if(nitmax.eq.3) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 0

      elseif(nitmax.eq.4) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 1
       iskip_lu(4, nq) = 0

      elseif(nitmax.eq.5) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 1
       iskip_lu(4, nq) = 2
       iskip_lu(5, nq) = 0

      elseif(nitmax.eq.6) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 3
       iskip_lu(4, nq) = 1
       iskip_lu(5, nq) = 2
       iskip_lu(6, nq) = 0

      elseif(nitmax.eq.7) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 3
       iskip_lu(4, nq) = 4
       iskip_lu(5, nq) = 1
       iskip_lu(6, nq) = 2
       iskip_lu(7, nq) = 0

      elseif(nitmax.eq.8) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 5
       iskip_lu(4, nq) = 3
       iskip_lu(5, nq) = 4 
       iskip_lu(6, nq) = 1
       iskip_lu(7, nq) = 2 
       iskip_lu(8, nq) = 0

      elseif(nitmax.eq.9) then

       iskip_lu(1, nq) = 0
       iskip_lu(2, nq) = 0
       iskip_lu(3, nq) = 5
       iskip_lu(4, nq) = 3
       iskip_lu(5, nq) = 4
       iskip_lu(6, nq) = 6
       iskip_lu(7, nq) = 1
       iskip_lu(8, nq) = 2
       iskip_lu(9, nq) = 0

      elseif(nitmax.eq.10) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  = 7
       iskip_lu(4, nq)  = 3
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  = 6
       iskip_lu(7, nq)  = 1
       iskip_lu(8, nq)  = 2
       iskip_lu(9, nq)  = 4
       iskip_lu(10, nq) = 0

      elseif(nitmax.eq.11) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  = 8
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 3
       iskip_lu(6, nq)  = 4
       iskip_lu(7, nq)  = 6
       iskip_lu(8, nq)  = 1
       iskip_lu(9, nq)  = 7
       iskip_lu(10, nq) = 2
       iskip_lu(11, nq) = 0

      elseif(nitmax.eq.12) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  = 9
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 8
       iskip_lu(6, nq)  = 3
       iskip_lu(7, nq)  = 4
       iskip_lu(8, nq)  = 6
       iskip_lu(9, nq)  = 1
       iskip_lu(10, nq) = 2
       iskip_lu(11, nq) = 7
       iskip_lu(12, nq) = 0

      elseif(nitmax.eq.13) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  = 9
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 3
       iskip_lu(6, nq)  = 8
       iskip_lu(7, nq)  = 4
       iskip_lu(8, nq)  =10
       iskip_lu(9, nq)  = 1
       iskip_lu(10, nq) = 6
       iskip_lu(11, nq) = 2
       iskip_lu(12, nq) = 7
       iskip_lu(13, nq) = 0

      elseif(nitmax.eq.14) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =11
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  = 3
       iskip_lu(7, nq)  = 8
       iskip_lu(8, nq)  = 4
       iskip_lu(9, nq)  = 6
       iskip_lu(10, nq) = 1
       iskip_lu(11, nq) =10
       iskip_lu(12, nq) = 2
       iskip_lu(13, nq) = 7
       iskip_lu(14, nq) = 0

      elseif(nitmax.eq.15) then

c       iskip_lu(nitmax, 1)  = 0
c       iskip_lu(nitmax, 2)  = 0
c       iskip_lu(nitmax, 3)  = 2 
c       iskip_lu(nitmax, 4)  = 4
c       iskip_lu(nitmax, 5)  = 6
c       iskip_lu(nitmax, 6)  = 8 
c       iskip_lu(nitmax, 7)  = 10 
c       iskip_lu(nitmax, 8)  =12
c       iskip_lu(nitmax, 9)  =11
c       iskip_lu(10, nq)  = 9
c       iskip_lu(11, nq)  = 7
c       iskip_lu(12, nq)  = 5 
c       iskip_lu(13, nq)  = 3
c       iskip_lu(14, nq)  =    1
c       iskip_lu(15, nq)  = 0

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =11
       iskip_lu(4, nq)  = 5 
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  = 3
       iskip_lu(7, nq)  = 8 
       iskip_lu(8, nq)  = 4
       iskip_lu(9, nq)  =12
       iskip_lu(10, nq) = 6
       iskip_lu(11, nq) = 1
       iskip_lu(12, nq) =10 
       iskip_lu(13, nq) = 2 
       iskip_lu(14, nq) = 7
       iskip_lu(15, nq) = 0

      elseif(nitmax.eq.16) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =11
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  = 3
       iskip_lu(7, nq)  = 8
       iskip_lu(8, nq)  =13
       iskip_lu(9, nq)  = 4
       iskip_lu(10, nq) =12
       iskip_lu(11, nq) = 6
       iskip_lu(12, nq) = 1
       iskip_lu(13, nq) =10
       iskip_lu(14, nq) = 2
       iskip_lu(15, nq) = 7
       iskip_lu(16, nq) = 0

      elseif(nitmax.eq.17) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =14
       iskip_lu(4, nq)  = 8
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  =13
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  = 9
       iskip_lu(10, nq) = 4
       iskip_lu(11, nq) =12
       iskip_lu(12, nq) = 6
       iskip_lu(13, nq) = 1
       iskip_lu(14, nq) =10
       iskip_lu(15, nq) = 2
       iskip_lu(16, nq) = 7
       iskip_lu(17, nq) = 0

      elseif(nitmax.eq.18) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =15
       iskip_lu(4, nq)  =10
       iskip_lu(5, nq)  = 3
       iskip_lu(6, nq)  =13
       iskip_lu(7, nq)  = 5
       iskip_lu(8, nq)  = 9
       iskip_lu(9, nq)  =12
       iskip_lu(10, nq) = 6
       iskip_lu(11, nq) = 4
       iskip_lu(12, nq) = 8
       iskip_lu(13, nq) = 1
       iskip_lu(14, nq) =11
       iskip_lu(15, nq) = 2
       iskip_lu(16, nq) =14
       iskip_lu(17, nq) = 7
       iskip_lu(18, nq) = 0

      elseif(nitmax.eq.19) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =14
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  =13
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  = 8
       iskip_lu(10, nq) =16
       iskip_lu(11, nq) = 4
       iskip_lu(12, nq) =12
       iskip_lu(13, nq) = 6
       iskip_lu(14, nq) = 1
       iskip_lu(15, nq) =10
       iskip_lu(16, nq) = 2
       iskip_lu(17, nq) =15
       iskip_lu(18, nq) = 7
       iskip_lu(19, nq) = 0

      elseif(nitmax.eq.20) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =14
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  =16
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  = 8
       iskip_lu(10, nq) =17
       iskip_lu(11, nq) = 4
       iskip_lu(12, nq) =12
       iskip_lu(13, nq) = 6
       iskip_lu(14, nq) = 1
       iskip_lu(15, nq) =10
       iskip_lu(16, nq) =15
       iskip_lu(17, nq) = 2
       iskip_lu(18, nq) = 7
       iskip_lu(19, nq) =13
       iskip_lu(20, nq) = 0

      elseif(nitmax.eq.21) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =14
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  =16
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  =12
       iskip_lu(10, nq) = 8
       iskip_lu(11, nq) =18
       iskip_lu(12, nq) = 4
       iskip_lu(13, nq) =12
       iskip_lu(14, nq) = 6
       iskip_lu(15, nq) = 1
       iskip_lu(16, nq) =10
       iskip_lu(17, nq) =17
       iskip_lu(18, nq) = 2
       iskip_lu(19, nq) = 7
       iskip_lu(20, nq) =15
       iskip_lu(21, nq) = 0

      elseif(nitmax.eq.22) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =18
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  =16
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =14
       iskip_lu(9, nq)  = 8
       iskip_lu(10, nq) =13
       iskip_lu(11, nq) = 4
       iskip_lu(12, nq) =19
       iskip_lu(13, nq) = 6
       iskip_lu(14, nq) =12
       iskip_lu(15, nq) = 1
       iskip_lu(16, nq) =10
       iskip_lu(17, nq) =17
       iskip_lu(18, nq) = 2
       iskip_lu(19, nq) = 7
       iskip_lu(20, nq) =15
       iskip_lu(21, nq) =11
       iskip_lu(22, nq) = 0

      elseif(nitmax.eq.23) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =18
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  = 9
       iskip_lu(6, nq)  =16
       iskip_lu(7, nq)  = 3
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  =14
       iskip_lu(10, nq) = 8
       iskip_lu(11, nq) =13
       iskip_lu(12, nq) = 4
       iskip_lu(13, nq) =20
       iskip_lu(14, nq) =12
       iskip_lu(15, nq) = 6
       iskip_lu(16, nq) = 1
       iskip_lu(17, nq) =19
       iskip_lu(18, nq) =10
       iskip_lu(19, nq) =17
       iskip_lu(20, nq) = 2
       iskip_lu(21, nq) = 7
       iskip_lu(22, nq) =15
       iskip_lu(23, nq) = 0

      elseif(nitmax.eq.24) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =20
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  =14
       iskip_lu(6, nq)  = 9
       iskip_lu(7, nq)  =16
       iskip_lu(8, nq)  = 3
       iskip_lu(9, nq)  =11
       iskip_lu(10, nq) =18
       iskip_lu(11, nq) = 8
       iskip_lu(12, nq) =13
       iskip_lu(13, nq) = 4
       iskip_lu(14, nq) =21
       iskip_lu(15, nq) =12
       iskip_lu(16, nq) =19
       iskip_lu(17, nq) = 1
       iskip_lu(18, nq) =10
       iskip_lu(19, nq) = 6
       iskip_lu(20, nq) =17
       iskip_lu(21, nq) = 2
       iskip_lu(22, nq) = 7
       iskip_lu(23, nq) =15
       iskip_lu(24, nq) = 0

      elseif(nitmax.eq.25) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =20
       iskip_lu(4, nq)  = 5
       iskip_lu(5, nq)  =14
       iskip_lu(6, nq)  = 9
       iskip_lu(7, nq)  =16
       iskip_lu(8, nq)  =11
       iskip_lu(9, nq)  =18
       iskip_lu(10, nq) = 8
       iskip_lu(11, nq) =13
       iskip_lu(12, nq) =22
       iskip_lu(13, nq) = 4
       iskip_lu(14, nq) =12
       iskip_lu(15, nq) = 6
       iskip_lu(16, nq) =19
       iskip_lu(17, nq) = 1
       iskip_lu(18, nq) =10
       iskip_lu(19, nq) =17
       iskip_lu(20, nq) = 2
       iskip_lu(21, nq) =20
       iskip_lu(22, nq) = 7
       iskip_lu(23, nq) =14
       iskip_lu(24, nq) = 3
       iskip_lu(25, nq) = 0

      elseif(nitmax.eq.26) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =21
       iskip_lu(4, nq)  =14
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  = 9
       iskip_lu(7, nq)  =16
       iskip_lu(8, nq)  = 3
       iskip_lu(9, nq)  =11
       iskip_lu(10, nq) =18
       iskip_lu(11, nq) = 8
       iskip_lu(12, nq) =13
       iskip_lu(13, nq) =23
       iskip_lu(14, nq) =12
       iskip_lu(15, nq) =20
       iskip_lu(16, nq) = 6
       iskip_lu(17, nq) =19
       iskip_lu(18, nq) = 1
       iskip_lu(19, nq) =10
       iskip_lu(20, nq) =17
       iskip_lu(21, nq) = 2
       iskip_lu(22, nq) =22
       iskip_lu(23, nq) = 4
       iskip_lu(24, nq) = 7
       iskip_lu(25, nq) =15
       iskip_lu(26, nq) = 0

      elseif(nitmax.eq.27) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =23
       iskip_lu(4, nq)  =14
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  = 9
       iskip_lu(7, nq)  =16
       iskip_lu(8, nq)  = 3
       iskip_lu(9, nq)  =11
       iskip_lu(10, nq) =18
       iskip_lu(11, nq) = 8
       iskip_lu(12, nq) =13
       iskip_lu(13, nq) =24
       iskip_lu(14, nq) = 4
       iskip_lu(15, nq) =20
       iskip_lu(16, nq) =22
       iskip_lu(17, nq) = 6
       iskip_lu(18, nq) =19
       iskip_lu(19, nq) = 1
       iskip_lu(20, nq) =10
       iskip_lu(21, nq) =17
       iskip_lu(22, nq) = 2
       iskip_lu(23, nq) =12
       iskip_lu(24, nq) =21
       iskip_lu(25, nq) = 7
       iskip_lu(26, nq) =15
       iskip_lu(27, nq) = 0

      elseif(nitmax.eq.28) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =24
       iskip_lu(4, nq)  =14
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  =21
       iskip_lu(7, nq)  = 9
       iskip_lu(8, nq)  =16
       iskip_lu(9, nq)  =11
       iskip_lu(10, nq) =18
       iskip_lu(11, nq) =20
       iskip_lu(12, nq) = 8
       iskip_lu(13, nq) = 4
       iskip_lu(14, nq) =25
       iskip_lu(15, nq) =12
       iskip_lu(16, nq) =23
       iskip_lu(17, nq) = 6
       iskip_lu(18, nq) =19
       iskip_lu(19, nq) = 1
       iskip_lu(20, nq) =10
       iskip_lu(21, nq) =17
       iskip_lu(22, nq) = 2
       iskip_lu(23, nq) =22
       iskip_lu(24, nq) = 3
       iskip_lu(25, nq) =13
       iskip_lu(26, nq) = 7
       iskip_lu(27, nq) =15
       iskip_lu(28, nq) = 0

      elseif(nitmax.eq.29) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =25
       iskip_lu(4, nq)  =14
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  =21
       iskip_lu(7, nq)  = 9
       iskip_lu(8, nq)  =16
       iskip_lu(9, nq)  = 3
       iskip_lu(10, nq) =11
       iskip_lu(11, nq) =18
       iskip_lu(12, nq) = 8
       iskip_lu(13, nq) =13
       iskip_lu(14, nq) =26
       iskip_lu(15, nq) =20
       iskip_lu(16, nq) =12
       iskip_lu(17, nq) =23
       iskip_lu(18, nq) = 6
       iskip_lu(19, nq) =19
       iskip_lu(20, nq) = 1
       iskip_lu(21, nq) = 4 
       iskip_lu(22, nq) =24
       iskip_lu(23, nq) =10
       iskip_lu(24, nq) =17
       iskip_lu(25, nq) = 2
       iskip_lu(26, nq) =22
       iskip_lu(27, nq) = 7
       iskip_lu(28, nq) =15
       iskip_lu(29, nq) = 0

      elseif(nitmax.eq.30) then

       iskip_lu(1, nq)  = 0
       iskip_lu(2, nq)  = 0
       iskip_lu(3, nq)  =25
       iskip_lu(4, nq)  =14
       iskip_lu(5, nq)  = 5
       iskip_lu(6, nq)  =21
       iskip_lu(7, nq)  = 9
       iskip_lu(8, nq)  =16
       iskip_lu(9, nq)  = 3
       iskip_lu(10, nq) =11
       iskip_lu(11, nq) =18
       iskip_lu(12, nq) =24
       iskip_lu(13, nq) = 8
       iskip_lu(14, nq) =13 
       iskip_lu(15, nq) =27
       iskip_lu(16, nq) = 4
       iskip_lu(17, nq) =20
       iskip_lu(18, nq) =12 
       iskip_lu(19, nq) = 6
       iskip_lu(20, nq) =19
       iskip_lu(21, nq) = 1
       iskip_lu(22, nq) =17
       iskip_lu(23, nq) =23
       iskip_lu(24, nq) =10
       iskip_lu(25, nq) =26
       iskip_lu(26, nq) = 2
       iskip_lu(27, nq) =22
       iskip_lu(28, nq) = 7
       iskip_lu(29, nq) =15
       iskip_lu(30, nq) = 0

      else
       iskip_lu(:,:) = 0
      endif

      end
