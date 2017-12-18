c***********************************************************************
c     $Date: 2011-06-20 11:10:00 +0200 (lun 20 jun 2011) $
c     $Revision: 40 $
c     $Author: NicolasAlferez $
c***********************************************************************
      subroutine move( rot, param_real, translation_pur)
c***********************************************************************
c_P                          O N E R A
c
c_DC  DATE_C : Janvier-2011  --  ALFEREZ
c
c_U   USER : ALFEREZ
c
c     ACT
c_A   Appel de la routine de calcul du torseur du mouvement relatif
c_A   du maillage
c
c
c     VAL
c
c     INP
c_I    temps  :  temps a l'iteration n (pas de temps)
c_I              ou n+1/2 
c
c     OUT
c_O    /relative/
c_O    rot    :  vecteur rotation du maillage
c_O    omega  :  vecteur vitesse angulaire du maillage
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      REAL_E  param_real(0:*), rot(4,3)
      INTEGER_E translation_pur

C Var loc
      REAL_E axe(3),axenorm,tnx,tny,tnz,teta,tetap,
     & racperpair11,racperpair21,racperpair31,racperpair12,
     & racperpair22,racperpair32,racperpair13,racperpair23,racperpair33,
     & racperimp11,racperimp21,racperimp31,racperimp12,
     & racperimp22,racperimp32,racperimp13,racperimp23,racperimp33,
     & teta_m1,teta_p1,f,fp,t,pi, teta_out

c.....formule de mvt
c      f(t)  = param_real(ROT_AMPLI)*sin(param_real(ROT_FREQ)*t)
c      fp(t) = param_real(ROT_AMPLI)*param_real(ROT_FREQ)
c     &        *cos(param_real(ROT_FREQ)*t)

c      f(t)      = param_real(ROT_AMPLI)*(1.-cos(param_real(ROT_FREQ)*t))
c      fp(t)     = param_real(ROT_AMPLI)*param_real(ROT_FREQ)
c     &           *sin(param_real(ROT_FREQ)*t)

      !pi    = 4.*atan(1.)

      axe(1)=param_real(ROT_OMEGA)
      axe(2)=param_real(ROT_OMEGA+1)
      axe(3)=param_real(ROT_OMEGA+2)

      translation_pur=1

      axenorm = sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
      if(axenorm.le.1.e-7)then
         tnx  =0.
         tny  =0.
         tnz  =0.
         teta =0.
         tetap=0.
         translation_pur=0
         goto 111
      endif
 
      tnx = axe(1) / axenorm
      tny = axe(2) / axenorm
      tnz = axe(3) / axenorm

      !calcul de la duree par rapport a l'instant de mise en route de l'ale
      !tps      = tt - param_real(ROT_INIT)
 
      !teta_out = f(tps)
      !tetap    = fp(tps)
      teta_out = param_real(ROT_TETA)
      tetap    = param_real(ROT_TETAP)
      teta     = teta_out
      !teta     = teta_out - teta_init_ale

c       write(*,'(a,4f17.13)')'teta,tetap,teta-teta_init ', 
c     & param_real(ROT_AMPLI)*param_real(ROT_FREQ), 
c     & cos(param_real(ROT_FREQ)*tps), param_real(ROT_FREQ)*tps, tps

      !pi    = 4.*atan(1.)

c       if(iverbs.ge.1.and.mpiproc.le.1.and.ldisplay)
c       write(*,'(a,4f17.8)')'teta,tetap,teta-teta_init ', 
c     & (teta_out)/pi*180., tetap/pi*180., teta/pi*180.,
c     & tps*param_real(ROT_FREQ)/2./pi

c######## Matrice de rotation d'angle teta 
c
  111 continue

       racperpair11 = tnx*tnx * (1.-cos(teta)) + cos(teta)
       racperpair21 = tnx*tny * (1.-cos(teta))
       racperpair31 = tnx*tnz * (1.-cos(teta))
       racperpair12 = tny*tnx * (1.-cos(teta))
       racperpair22 = tny*tny * (1.-cos(teta)) + cos(teta)
       racperpair32 = tny*tnz * (1.-cos(teta))
       racperpair13 = tnz*tnx * (1.-cos(teta))
       racperpair23 = tnz*tny * (1.-cos(teta))
       racperpair33 = tnz*tnz * (1.-cos(teta)) + cos(teta)
 
       racperimp11  =  0.
       racperimp21  =  tnz * sin(teta)
       racperimp31  = -tny * sin(teta)
       racperimp12  = -tnz * sin(teta)
       racperimp22  =  0.
       racperimp32  =  tnx * sin(teta)
       racperimp13  =  tny * sin(teta)
       racperimp23  = -tnx * sin(teta)
       racperimp33  =  0.

c########

       rot(1,1)  = racperpair11 + racperimp11
       rot(2,1)  = racperpair21 + racperimp21
       rot(3,1)  = racperpair31 + racperimp31
       rot(1,2)  = racperpair12 + racperimp12
       rot(2,2)  = racperpair22 + racperimp22
       rot(3,2)  = racperpair32 + racperimp32
       rot(1,3)  = racperpair13 + racperimp13
       rot(2,3)  = racperpair23 + racperimp23
       rot(3,3)  = racperpair33 + racperimp33
       

       rot(4,1) = tnx*tetap
       rot(4,2) = tny*tetap
       rot(4,3) = tnz*tetap
      
c      teta = teta*180./(4.*atan(1.))
c      open(545,file='teta',form='formatted',status='old')

c       tetatot=tetatot+teta
c       tetadebug=tetatot
c       write(545,1)temps,tetatot*180/3.14159256,teta*180/3.1415925,
c     &             nitrun,lteta
c    &       omegar(1),omegar(2),omegar(3),tt 
      
c      write(*,'(a,3f12.7)')'x ',rot(1,1),rot(1,2),rot(1,3)
c      write(*,'(a,3f12.7)')'y ',rot(2,1),rot(2,2),rot(2,3)
c      write(*,'(a,3f12.7)')'z ',rot(3,1),rot(3,2),rot(3,3)

c 1     format (3(e15.8,2x),I3,2x,L5)

       end
