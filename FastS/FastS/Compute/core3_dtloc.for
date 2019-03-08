c***********************************************************************
c     $Date: 2013-08-26 16:00:23 +0200 (lun. 26 août 2013) $
c     $Revision: 58 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine core3_dtloc(ndom, nitcfg,  param_int, param_real,
     &                     ind_loop,
     &                     rop_1, rop, drodm, coe)
c***********************************************************************
c_U   USER : PECHIER
c
c     ACT
c_A    Mise a jour du residu explicite
c
c     VAL
c
c***********************************************************************
      implicit none

#include "FastS/param_solver.h"

      REAL_E cgamma1, cgamma2, cgamma3,cgamma4
      parameter( cgamma1 =  0.5             )
      parameter( cgamma2 =  0.9106836025229591 )
      parameter( cgamma3 =  0.3660254037844387 )
      parameter( cgamma4 =  3.6427344100918355 )

      include 'omp_lib.h'

      INTEGER_E ndom,nitcfg,ind_loop(6), param_int(0:*)
 
      REAL_E     rop( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   rop_1( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E   drodm( param_int(NDIMDX) * param_int(NEQ) )
      REAL_E     coe( param_int(NDIMDX) * param_int(NEQ_COE) )

      REAL_E param_real(0:*)
C Var loc
      INTEGER_E incmax,l,i,j,k,lij,v1,v2,v3,v4,v5,v6,b,ind
      REAL_E coefW,ro_old,u_old,v_old,w_old,t_old,nu_old,r_1,cvinv,
     &  roe_old, cv, cvinv2, tab1,tab2,tab3,tab4,tab5,tab6,
     & rho,rhou,rhov,rhoe,anulam,temp01,cmus1,coesut,t1
      REAL_E a1,a2,a3,a4,a5,a6,a7,a1_,a2_,a3_,a4_,a5_,a6_,a7_
      INTEGER_E b1,b2,b3,b4,b5,b6,b7,shift

#include "FastS/formule_param.h"

      
      cv     = param_real( CVINF)   ! Cv
      cvinv  = 1./cv   
      cvinv2 = 0.5*cvinv

      v1 = 0
      v2 =   param_int(NDIMDX)
      v3 = 2*param_int(NDIMDX)
      v4 = 3*param_int(NDIMDX)
      v5 = 4*param_int(NDIMDX)
      v6 = 5*param_int(NDIMDX)

      b=param_int(RK)

      ind=param_int(NSSITER)/param_int(LEVEL)
      !shift=param_int(SHIFTLOCAL)
      shift=771840
      !print*, param_int(NDIMDX)

     

c---------------------------------------------------------------------
c-----Determination des coefficients : Explicit global
c---------------------------------------------------------------------
      if (param_int(EXPLOC)==0) then            

         if (b==1) then ! Runge-Kutta ordre 1
            if (nitcfg.eq.1) then
               coefW = 1.
            end if
         elseif (b==2) then ! Runge-Kutta ordre 2
            if (nitcfg.eq.1) then
               coefW = 1.
            elseif (nitcfg.eq.2) then
               coefW = 0.5
            end if
         elseif (b==3) then ! Runge-Kutta ordre 3
            if (nitcfg.eq.1) then
               coefW = cgamma1
            elseif (nitcfg.eq.2) then
               coefW = cgamma2
            elseif (nitcfg.eq.3) then
               coefW = cgamma3
            end if
         elseif (b==5) then ! Runge-Kutta 5 ordre 4
            if (nitcfg.eq.1) then
               !coefW = 1.3512072563171387
               !coefW = 9.7618354692056E-2
               coefW = 0.1496590219993
            elseif (nitcfg.eq.2) then
               !coefW = -9.7900167107582E-002
               !coefW = 0.4122532929155
               coefW = 0.3792103129999
            elseif (nitcfg.eq.3) then
               !coefW = -1.7024143934249878
               !coefW = 0.4402169639311
               coefW = 0.8229550293869
            elseif (nitcfg.eq.4) then
               !coefW = -0.35120719671249390
               !coefW = 1.426311463224
               coefW = 0.6994504559488
            elseif (nitcfg.eq.5) then
               !coefW = 0.67560362815856934
               !coefW = 0.1978760537318
               coefW = 0.1530572479681
            end if
         elseif (b==12) then ! Runge-Kutta 5 ordre 4
            if (nitcfg.eq.1) then
               !coefW = 0.0367762454319673
               coefW = 0.0650008435125904
            elseif (nitcfg.eq.2) then
               !coefW = 0.3136296607553959
               coefW = 0.0161459902249842
            elseif (nitcfg.eq.3) then
               !coefW = 0.1531848691869027
               coefW = 0.5758627178358159
            elseif (nitcfg.eq.4) then
               !coefW = 0.0030097086818182
               coefW = 0.1649758848361671
            elseif (nitcfg.eq.5) then
               !coefW = 0.3326293790646110
               coefW = 0.3934619494248182
            elseif (nitcfg.eq.6) then
               !coefW = 0.2440251405350864
               coefW = 0.0443509641602719
            elseif (nitcfg.eq.7) then
               !coefW = 0.3718879239592277 
               coefW = 0.2074504268408778
            elseif (nitcfg.eq.8) then
               !coefW = 0.6204126221582444
               coefW = 0.6914247433015102
            elseif (nitcfg.eq.9) then
               !coefW = 0.1524043173028741
               coefW = 0.3766646883450449 
            elseif (nitcfg.eq.10) then
               !coefW = 0.0760894927419266
               coefW = 0.0757190350155483  
            elseif (nitcfg.eq.11) then
               !coefW = 0.0077604214040978
               coefW = 0.2027862031054088
            elseif (nitcfg.eq.12) then
               !coefW = 0.0024647284755382
               coefW = 0.2167029365631842
            elseif (nitcfg.eq.13) then
               !coefW = 0.0780348340049386
            elseif (nitcfg.eq.14) then
               !coefW = 5.5059777270269628
            end if
         elseif (b==4) then ! Runge-Kutta ordre 4
            if (nitcfg.eq.1) then
               a1=0.5
               a2=0.
               a3=0.
               a4=0.
               b1=0
               b2=0
               b3=0
               b4=0
            elseif (nitcfg.eq.2) then
               a1=-0.5
               a2=0.5
               a3=0.
               a4=0.
               b1=-1
               b2=0
               b3=0
               b4=0
            elseif (nitcfg.eq.3) then
               a1=-0.5
               a2=1.0
               a3=0.
               a4=0.
               b1=-1
               b2=0
               b3=0
               b4=0
            elseif (nitcfg.eq.4) then
               a1=1./6.
               a2=1./3.
               a3=-2./3.
               a4=1./6.
               b1=-3
               b2=-2
               b3=-1
               b4=0
            end if

c$$$         elseif (b==5) then ! Runge-Kutta ordre 6
c$$$            if (nitcfg.eq.1) then
c$$$               a1=0.37726891533136803 !1./6.
c$$$               a2=0.!0.
c$$$               a3=0.!0.
c$$$               a4=0.!0.
c$$$               a5=0.!0.
c$$$               a6=0.!0.
c$$$               a7=0.!0.
c$$$               b1=0!0!0!0!0
c$$$               b2=0!0!0!0!0
c$$$               b3=0!0!0!0!0
c$$$               b4=0!0!0!0!0
c$$$               b5=0!0
c$$$               b6=0!0
c$$$               b7=0!0
c$$$            elseif (nitcfg.eq.2) then
c$$$               a1=0.!-0.38294237852096558 
c$$$               a2=0.37726891533136803! 0.54960906505584717
c$$$               a3=0.!0.
c$$$               a4=0.!0.
c$$$               a5=0.!0.
c$$$               a6=0.!0.
c$$$               a7=0.!0.
c$$$               b1=-1!0!0!0
c$$$               b2=0!0!0!0
c$$$               b3=0!0!0!0
c$$$               b4=0!0!0!0
c$$$               b5=0
c$$$               b6=0
c$$$               b7=0 
c$$$            elseif (nitcfg.eq.3) then
c$$$               a1=-0.13427369479397216! 0.30110452324151993
c$$$               a2=-0.13427369479397216!-0.50798253342509270
c$$$               a3= 0.24299522053739600 ! 0.37354466319084167 
c$$$               a4=0.!0.
c$$$               a5=0.!0.
c$$$               a6=0.!0.
c$$$               a7=0.!0.
c$$$               b1=-2!-2
c$$$               b2=-1!-1
c$$$               b3=0!0
c$$$               b4=0!0
c$$$               b5=0!0
c$$$               b6=0!0
c$$$               b7=0!0  
c$$$            elseif (nitcfg.eq.4) then
c$$$               a1=-8.9406152842269354E-002  !-0.17133980244398117
c$$$               a2=-8.9406152842269354E-002!0.33792908117175102
c$$$               a3=-8.9406152842269410E-002!-0.35600895434617996 
c$$$               a4= 0.23845893284628999   !0.35608631372451782 
c$$$               a5=0.!0.
c$$$               a6=0.!0.
c$$$               a7=0.!0.
c$$$               b1=-3!-3
c$$$               b2=-2!-2
c$$$               b3=-1!-1
c$$$               b4=0!0
c$$$               b5=0!0
c$$$               b6=0!0
c$$$               b7=0!0  
c$$$            elseif (nitcfg.eq.5) then
c$$$               a1= 5.3144953169677961E-002!-3.8486570119857788E-002
c$$$               a2= 5.3144953169677961E-002!0.34739521145820618
c$$$               a3=-3.6491815853282469E-002!-0.40116741508245468
c$$$               a4=-5.6656372726150567E-002!-6.1162561178207397E-002
c$$$               a5=0.28763214630840800 !0.32008802890777588
c$$$               a6=0.
c$$$               a7=0.
c$$$               b1=-4
c$$$               b2=-3
c$$$               b3=-2
c$$$               b4=-1
c$$$               b5=0
c$$$               b6=0
c$$$               b7=0  
c$$$            elseif (nitcfg.eq.6) then
c$$$               a1=0.20392320305109024
c$$$               a2=-0.57157905399799347  
c$$$               a3=0.47288817912340164 
c$$$               a4=0.21581980586051941 
c$$$               a5=-0.62545979022979736
c$$$               a6=0.47107434272766113
c$$$               a7=0.
c$$$               b1=-5
c$$$               b2=-4
c$$$               b3=-3
c$$$               b4=-2
c$$$               b5=-1
c$$$               b6=0
c$$$               b7=0  
c$$$            elseif (nitcfg.eq.7) then
c$$$               a1=-596./5280.
c$$$               a2=9./11.
c$$$               a3=-1332./1760.
c$$$               a4=-423./440.
c$$$               a5=-4./15.
c$$$               a6=196./165.
c$$$               a7=11./120.
c$$$               b1=-6
c$$$               b2=-5
c$$$               b3=-4
c$$$               b4=-3
c$$$               b5=-2
c$$$               b6=-1
c$$$               b7=0  
 
           ! end if

        !print*, a1+a2+a3+a4+a5+a6+a7
         
        end if
c---------------------------------------------------------------------
c-----Determination des coefficients : Explicit local (RK 2)
c---------------------------------------------------------------------

      elseif (param_int(EXPLOC)==1.and.param_int(RK)==2) then 

       if (ind==2) then ! zones de + gd niveau en temps 

         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
            !0.25 pour le slow buffer
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
         end if

       else ! zones dont le niveau en temps n'est pas le + gd
         
         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,'coefW= ', coefW
         elseif (MOD(nitcfg,ind)==ind/2) then 
            !coefW = 0.25/float(param_int(LEVEL))
            coefW = 0
            !print*, 'coefW= ', coefW
         end if


       end if



c---------------------------------------------------------------------
c-----Determination des coefficients : Tang & Warnecke
c---------------------------------------------------------------------

      elseif (param_int(EXPLOC)==5.and.param_int(RK)==2) then 

       if (ind==2) then ! zones de + gd niveau en temps 

         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
            !0.25 pour le slow buffer
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
         end if

       else ! zones dont le niveau en temps n'est pas le + gd
         
         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, 'coefW= ', coefW
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,'coefW= ', coefW
         elseif (MOD(nitcfg,ind)==ind/2) then 
            coefW = 0.25/float(param_int(LEVEL))
            !coefW = 0
            !print*, 'coefW= ', coefW
         end if


       end if


c-----------------------Constantinescu----------------------------------


      elseif (param_int(EXPLOC)==2.and.param_int(RK)==2) then 
         
       if (ind==2) then ! zones de + gd niveau en temps 

         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
            !0.25 pour le slow buffer
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,coefW
         end if

       else ! zones dont le niveau en temps n'est pas le + gd

         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
            !0.25 pour le slow buffer
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,coefW
         elseif (MOD(nitcfg,ind)==ind/2) then 
            coefW = 1./float(param_int(LEVEL))
            ! print*, coefW
         elseif (MOD(nitcfg,ind)==ind/2+1) then 
            coefW = 1./float(param_int(LEVEL))
            ! print*, coefW
         end if

       end if



c--------------------Constantinescu ac rk3----------------------------


      else if (param_int(EXPLOC)==3.and.param_int(RK)==2) then


         if (ind==3) then ! zones de + gd niveau en temps

            if (MOD(nitcfg,ind).eq.0) then
              coefW = cgamma3/float(param_int(LEVEL))
              !print*, 'coefW= ',coefW
            elseif (MOD(nitcfg,ind).eq.1) then
              coefW = cgamma1/float(param_int(LEVEL))
              !print*, 'coefW= ', coefW
            elseif (MOD(nitcfg,ind).eq.2) then
              coefW = cgamma2/float(param_int(LEVEL))
              !print*, 'coefW= ', coefW
            end if
            
       
         else  !zones dont le niveau en temps n'est pas le + gd 

            if (MOD(nitcfg,ind).eq.0) then
              coefW = cgamma3/float(param_int(LEVEL))    
            elseif (MOD(nitcfg,ind).eq.1) then
              coefW = cgamma1/float(param_int(LEVEL))                           
            elseif (MOD(nitcfg,ind).eq.2) then
              coefW = cgamma2/float(param_int(LEVEL))  
            elseif (MOD(nitcfg,ind).eq.ind/2) then
              coefW = cgamma3/float(param_int(LEVEL))
            elseif (MOD(nitcfg,ind).eq.ind/2+1) then
              coefW = cgamma1/float(param_int(LEVEL)) 
            elseif (MOD(nitcfg,ind).eq.ind/2+2) then
              coefW = cgamma2/float(param_int(LEVEL)) 
            end if

         end if  



c-------------------- RK3 local------------------------------------
      elseif (param_int(EXPLOC)==2.and.param_int(RK)==3) then 
         
         if (MOD(nitcfg,ind)==1) then
            coefW = cgamma1/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         elseif (MOD(nitcfg,ind)==ind/2) then
            coefW = cgamma2/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         elseif (MOD(nitcfg,ind)==ind-1) then
            coefW = cgamma3/float(param_int(LEVEL))
            !print*,'coeff, ndom = ',coefW, ndom
         else if (MOD(nitcfg,ind)==ind/4.and.ind.ne.4) then
            coefW = cgamma4/float(param_int(LEVEL))
            !print*,'coeff= ',coefW
                    
         end if 


c--------------------RK2 local test--------------------------------

      elseif (param_int(EXPLOC)==4.and.param_int(RK)==2) then 

       if (ind==2) then ! zones de + gd niveau en temps 

         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
            !0.25 pour le slow buffer
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,coefW
         end if

       else ! zones dont le niveau en temps n'est pas le + gd
         
         if (MOD(nitcfg,ind)==0) then               
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
         elseif (MOD(nitcfg,ind)==1) then 
            coefW = 1./float(param_int(LEVEL))
            !print*,coefW
         elseif (MOD(nitcfg,ind)==ind/2) then 
            coefW = 0.25/float(param_int(LEVEL))
            !print*, coefW
         elseif (MOD(nitcfg,ind)==ind-1) then 
            coefW = 0.5/float(param_int(LEVEL))
            !print*, coefW
         end if


       end if


   


      end if
c---------------------------------------------------------------------


      IF(param_int(NEQ).eq.5) THEN


         !print*, 'param(NDIMDX)= ', param_int(NDIMDX)

         If (param_int(ITYPZONE).ne.3) then

            if (b==4) then !!! En RK4


#include  "FastS/Compute/loop_core_begin.for"

          a1_    = a1*coe(l +v1)
          a2_    = a2*coe(l +v1)
          a3_    = a3*coe(l +v1)
          a4_    = a4*coe(l +v1)

          tab1 = drodm(l + b1*shift +v1)*a1_+drodm(l + b2*shift+v1)*a2_+
     & drodm(l + b3*shift + v1)*a3_+drodm(l + b4*shift + v1)*a4_  
          tab2 = drodm(l + b1*shift +v2)*a1_+drodm(l + b2*shift+v2)*a2_+
     & drodm(l + b3*shift + v2)*a3_+drodm(l + b4*shift + v2)*a4_  
          tab3 = drodm(l + b1*shift +v3)*a1_+drodm(l + b2*shift+v3)*a2_+
     & drodm(l + b3*shift + v3)*a3_+drodm(l + b4*shift + v3)*a4_  
          tab4 = drodm(l + b1*shift +v4)*a1_+drodm(l + b2*shift+v4)*a2_+
     & drodm(l + b3*shift + v4)*a3_+drodm(l + b4*shift + v4)*a4_  
          tab5 = drodm(l + b1*shift +v5)*a1_+drodm(l + b2*shift+v5)*a2_+
     & drodm(l + b3*shift + v5)*a3_+drodm(l + b4*shift + v5)*a4_  

c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          w_old       = rop(l +v4)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))

#include  "FastS/Compute/loop_end.for"

          else !!! Pas en RK4
               

#include  "FastS/Compute/loop_core_begin.for"

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!! Schema de Constantinescu !!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(param_int(EXPLOC)==2.and.param_int(RK)==2.and.
     &  MOD(nitcfg,ind).eq.0.and.ind.ne.b) then
 
        if( param_int(LEVEL).lt.param_int(LEVELD)
     &   .and.param_int(LEVELG).lt.param_int(LEVEL) ) then

           if(l.ge.lij+ind_loop(2)-ind_loop(1)-3) then
       
            coefW = 0.25/float(param_int(LEVEL))
           else  

            coefW = 0.5/float(param_int(LEVEL))
           end if
        end if

        if( param_int(LEVELG).ge.param_int(LEVEL)
     &    .and.param_int(LEVEL).ge.param_int(LEVELD) ) then

           if(l.le.lij+3) then          

           coefW = 0.25/float(param_int(LEVEL))
           else 

           coefW = 0.5/float(param_int(LEVEL))
           end if
        end if

      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab4 =  drodm(l +v4)*r_1
          tab5 =  drodm(l +v5)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          w_old       = rop(l +v4)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))

            ! if (k==100.and.j==100.and.l-lij==50) then
            !    print*, 'coeff=',coefW, ndom
            ! end if

          !if (ndom==53.and.nitcfg==3) then
          !   print*, rop_1(l+v5), k
          !end if

          !if () then
          !   print*, rop_1(l+v4)
          !end if

#include  "FastS/Compute/loop_end.for"
         
          end if !!! Condition RK4

         Else

            if (b==4) then !!! En RK4

#include  "FastS/Compute/loop_core_begin.for"

          a1_    = a1*coe(l +v1)
          a2_    = a2*coe(l +v1)
          a3_    = a3*coe(l +v1)
          a4_    = a4*coe(l +v1)

          tab1 = drodm(l + b1*shift +v1)*a1_+drodm(l + b2*shift+v1)*a2_+
     & drodm(l + b3*shift + v1)*a3_+drodm(l + b4*shift + v1)*a4_  
          tab2 = drodm(l + b1*shift +v2)*a1_+drodm(l + b2*shift+v2)*a2_+
     & drodm(l + b3*shift + v2)*a3_+drodm(l + b4*shift + v2)*a4_  
          tab3 = drodm(l + b1*shift +v3)*a1_+drodm(l + b2*shift+v3)*a2_+
     & drodm(l + b3*shift + v3)*a3_+drodm(l + b4*shift + v3)*a4_  
          tab5 = drodm(l + b1*shift +v5)*a1_+drodm(l + b2*shift+v5)*a2_+
     & drodm(l + b3*shift + v5)*a3_+drodm(l + b4*shift + v5)*a4_  

c

          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))

         !if (ndom==0 .and. l-lij==199.and.j==250) then
          !  print*, rop_1(l+v1)
          !  print*, drodm(l+b3*shift+v1),drodm(l+b4*shift+v1),b3,b4
         !end if

          !print*, 'coucou_rk4'

#include  "FastS/Compute/loop_end.for"


            else if (b==6) then !!! En RK6

#include  "FastS/Compute/loop_core_begin.for"

          a1_    = a1*coe(l +v1)
          a2_    = a2*coe(l +v1)
          a3_    = a3*coe(l +v1)
          a4_    = a4*coe(l +v1)
          a5_    = a5*coe(l +v1)
          a6_    = a6*coe(l +v1)
          a7_    = a7*coe(l +v1)

          tab1 = drodm(l + b1*shift +v1)*a1_+drodm(l + b2*shift+v1)*a2_+
     & drodm(l + b3*shift + v1)*a3_+drodm(l + b4*shift + v1)*a4_ +  
     & drodm(l + b5*shift +v1)*a5_+drodm(l + b6*shift+v1)*a6_+
     & drodm(l + b7*shift +v1)*a7_
 
          tab2 = drodm(l + b1*shift +v2)*a1_+drodm(l + b2*shift+v2)*a2_+
     & drodm(l + b3*shift + v2)*a3_+drodm(l + b4*shift + v2)*a4_  +
     & drodm(l + b5*shift +v2)*a5_+drodm(l + b6*shift+v2)*a6_+
     & drodm(l + b7*shift +v2)*a7_

          tab3 = drodm(l + b1*shift +v3)*a1_+drodm(l + b2*shift+v3)*a2_+
     & drodm(l + b3*shift + v3)*a3_+drodm(l + b4*shift + v3)*a4_  +  
     & drodm(l + b5*shift +v3)*a5_+drodm(l + b6*shift+v3)*a6_+
     & drodm(l + b7*shift +v3)*a7_
        
          tab5 = drodm(l + b1*shift +v5)*a1_+drodm(l + b2*shift+v5)*a2_+
     & drodm(l + b3*shift + v5)*a3_+drodm(l + b4*shift + v5)*a4_  +
     & drodm(l + b5*shift +v5)*a5_+drodm(l + b6*shift+v5)*a6_+
     & drodm(l + b7*shift +v5)*a7_

c

          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))


          !print*, 'coucou_rk4'
          !if (ndom==0 .and. l-lij==200.and.j==200) then
          !   print*, tab1, tab2, tab3, tab5
          !end if

#include  "FastS/Compute/loop_end.for"



          else !!!  Pas en RK4

#include  "FastS/Compute/loop_core_begin.for"


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!! Schema de Constantinescu !!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(param_int(EXPLOC)==2.and.param_int(RK)==2.and.
     &  MOD(nitcfg,ind).eq.0.and.ind.ne.b) then
 
        if( param_int(LEVEL).lt.param_int(LEVELD)
     &   .and.param_int(LEVELG).lt.param_int(LEVEL) ) then

           if(l.ge.lij+ind_loop(2)-ind_loop(1)-3) then
       
            coefW = 0.25/float(param_int(LEVEL))
           else  

            coefW = 0.5/float(param_int(LEVEL))
           end if
        end if

        if( param_int(LEVELG).ge.param_int(LEVEL)
     &    .and.param_int(LEVEL).ge.param_int(LEVELD) ) then

           if(l.le.lij+3) then          

           coefW = 0.25/float(param_int(LEVEL))
           else 

           coefW = 0.5/float(param_int(LEVEL))
           end if
        end if

      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !if (j==50.and.l.ge.lij+ind_loop(2)-ind_loop(1)-5) then
      !   print*, 'coefW= ', coefW," ",l-lij
      !end if

          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab5 =  drodm(l +v5)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.
 
          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))



          !if (ndom==0 .and. l-lij==100.and.j==200) then
          !   print*, tab1
          !end if


        !  print*, rop_1(l +v2)

      !if (ndom == 1.and.j==60.and.l-lij==86) then
      !   print*,'ro= ', drodm(l+v1), rop_1(l +v1)
      !end if
      !if (ndom == 0.and.j==59.and.l-lij==87) then !86
      !   print*,'ro= ',rop(l+v5), drodm(l +v5), rop_1(l +v5)
      !end if


#include  "FastS/Compute/loop_end.for"

          end if !!! Condition RK4

         Endif
      
      ELSE

      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

         If (param_int(ITYPZONE).ne.3) then

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab4 =  drodm(l +v4)*r_1
          tab5 =  drodm(l +v5)*r_1
          tab6 =  drodm(l +v6)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          w_old       = rop(l +v4)
          t_old       = rop(l +v5)
          nu_old      = rop(l +v6)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = (ro_old*w_old + tab4)*r_1


          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old
     &                                             +w_old*w_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3)
     &                              +rop_1(l +v4)*rop_1(l +v4))

          anulam=coesut*sqrt(rop_1(l+v5)*temp01)/(1.+cmus1/rop_1(l +v5))
          anulam=anulam*r_1*param_real(SA_REAL+SA_RATIOM-1) 

          t_old        = (ro_old*nu_old + tab6)*r_1
          t_old        = min(t_old,anulam)
          rop_1(l +v6) = max(t_old,0.)

#include  "FastS/Compute/loop_end.for"

         Else

#include  "FastS/Compute/loop_core_begin.for"
          r_1    = coefW*coe(l +v1)

          tab1 =  drodm(l +v1)*r_1
          tab2 =  drodm(l +v2)*r_1
          tab3 =  drodm(l +v3)*r_1
          tab5 =  drodm(l +v5)*r_1
          tab6 =  drodm(l +v6)*r_1
c
          ro_old      = rop(l +v1)
          u_old       = rop(l +v2)
          v_old       = rop(l +v3)
          t_old       = rop(l +v5)
          nu_old      = rop(l +v6)

          rop_1(l +v1)     = rop(l +v1) + tab1

          r_1            = 1./rop_1(l +v1)

          rop_1(l +v2)     = (ro_old*u_old + tab2)*r_1
          rop_1(l +v3)     = (ro_old*v_old + tab3)*r_1
          rop_1(l +v4)     = 0.

          roe_old        = ro_old*(cv*t_old + 0.5*( u_old*u_old
     &                                             +v_old*v_old))

          rop_1(l +v5)     = ( roe_old + tab5 )*r_1*cvinv
     &                    - cvinv2*( rop_1(l +v2)*rop_1(l +v2)
     &                              +rop_1(l +v3)*rop_1(l +v3))

          anulam=coesut*sqrt(rop_1(l+v5)*temp01)/(1.+cmus1/rop_1(l +v5))
          anulam=anulam*r_1*param_real(SA_REAL+SA_RATIOM-1) 

          t_old        = (ro_old*nu_old + tab6)*r_1
          t_old        = min(t_old,anulam)
          !rop_1(l +v6) = t_old
          rop_1(l +v6) = max(t_old,0.)

#include  "FastS/Compute/loop_end.for"


        Endif! typezone
      ENDIF! neq

      end
