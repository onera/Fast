c***********************************************************************
c     $Date: 2013-04-30 16:01:06 +0200 (mar. 30 avril 2013) $
c     $Revision: 56 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine init_rhs_dtloc(ndom, nitcfg, param_int, ndimdx, neq,
     &                   ind_loop,
     &                   drodm)
c***********************************************************************
c_P                          O N E R A
c
c_PR  PRODUIT :  mjdro3.sf/pechier1/e1/i1
c
c     ACT
c_A    Mise a jour de drodm apres calcul des flux.
c
c     INP
c_I    ndom   : numero-utilisateur du domaine
c_I    ipara  : direction du calcul des flux
c_I    neq    : nombre d equations du systeme
c_I    ndimdx : nombre de points maximal dans un domaine
c_I    flu    : flux aux interfaces
c
c     OUT
c
c     I/O
c_/    drodm  : increment des variables conservatives
c***********************************************************************
      implicit none

      REAL_E cxsi1, cxsi2,cxsi3,cxsi1_, cxsi2_,cxsi3_,cxsi4_,cxsi5_
      REAL_E cxsi6_,cxsi7_,cxsi8_,cxsi9_,cxsi10_,cxsi11_,cxsi12_,cxsi13_

      parameter( cxsi1   = -0.6830127018922193 )
      parameter( cxsi2   = - 4./3.             )
      parameter( cxsi3   = -0.96037658773652734 )

      !parameter( cxsi1_   = -1 )
      !parameter( cxsi2_   = -0.90871393680572510)
      !parameter( cxsi3_   = -4.8473219871520996)
      !parameter( cxsi4_   =  0.25992107391357422)

      !parameter( cxsi1_   = -0.4812317431372)
      !parameter( cxsi2_   = -1.049562606709)
      !parameter( cxsi3_   = -1.602529574275)
      !parameter( cxsi4_   = -1.778267193916)

      parameter( cxsi1_   = -0.4178904745)
      parameter( cxsi2_   = -1.192151694643)
      parameter( cxsi3_   = -1.697784692471)
      parameter( cxsi4_   = -1.514183444257)

      !parameter( cxsi1_   = -0.7188012108672410)
      !parameter( cxsi2_   = -0.7785331173421570)
      !parameter( cxsi3_   = -0.0053282796654044)
      !parameter( cxsi4_   = -0.8552979934029281)
      !parameter( cxsi5_   = -3.9564138245774565)
      !parameter( cxsi6_   = -1.5780575380587385)
      !parameter( cxsi7_   = -2.0837094552574054)
      !parameter( cxsi8_   = -0.7483334182761610)
      !parameter( cxsi9_   = -0.7032861106563359)
      !parameter( cxsi10_   = 0.0013917096117681)
      !parameter( cxsi11_   = -0.0932075369637460)
      !parameter( cxsi12_   = -0.9514200470875948)
      !parameter( cxsi13_   = -7.1151571693922548)

      !parameter( cxsi1_   = -0.0923311242368072)
      !parameter( cxsi2_   = -0.9441056581158819)
      !parameter( cxsi3_   = -4.3271273247576394)
      !parameter( cxsi4_   = -2.1557771329026072)
      !parameter( cxsi5_   = -0.9770727190189062)
      !parameter( cxsi6_   = -0.7581835342571139)
      !parameter( cxsi7_   = -1.7977525470825499)
      !parameter( cxsi8_   = -2.6915667972700770)
      !parameter( cxsi9_   = -4.6466798960268143)
      !parameter( cxsi10_  = -0.1539613783825189)
      !parameter( cxsi11_  = -0.5943293901830616)
      !parameter( cxsi12_   = -0.9514200470875948)
      !parameter( cxsi13_   = -7.1151571693922548)

#include "FastS/param_solver.h"

      INTEGER_E ndom, nitcfg, ndimdx, neq, param_int(0:*), ind_loop(6)

      REAL_E drodm( ndimdx * neq )
 

C Var loc
      INTEGER_E incmax, l, i,j,k,ne,lij,ltij,lt,b,ind, vg, lvo
      REAL_E ratio,coefH,xmut(1),rop(1) !!ajout pour feinter option de vecto

#include "FastS/formule_param.h"
#include "FastS/formule_mtr_param.h"


       if (param_int(EXPLOC)==1) then
         ind_loop(2)=ind_loop(2)+1
         ind_loop(1)=ind_loop(1)-1
      end if

      !print*, ind_loop(1) , ind_loop(2), ind_loop(3), ind_loop(4)


c$$$       if (param_int(EXPLOC)==2 .and. param_int(NIJK+4)==0) then
c$$$         ind_loop(2)=ind_loop(2)+1
c$$$         ind_loop(1)=ind_loop(1)-1
c$$$         ind_loop(4)=ind_loop(4)+1
c$$$         ind_loop(3)=ind_loop(3)-1
c$$$      else if(param_int(EXPLOC)==2 .and. param_int(NIJK+4).ne.0) then  
c$$$         ind_loop(2)=ind_loop(2)+1
c$$$         ind_loop(1)=ind_loop(1)-1
c$$$         ind_loop(4)=ind_loop(4)+1
c$$$         ind_loop(3)=ind_loop(3)-1
c$$$         ind_loop(6)=ind_loop(6)+1
c$$$         ind_loop(5)=ind_loop(5)-1         
c$$$      end if

      b   = param_int(RK)
      ind = param_int(NSSITER)/param_int(LEVEL)
     
c***********************************************************************
Calcul des coeffcients : Explicit global
c***********************************************************************

      if (param_int(EXPLOC)==0) then

      if (b.eq.3) then       ! Runge-Kutta ordre 3
       
         if (nitcfg.eq.2) then
             coefH = cxsi1
         else if (nitcfg.eq.3) then
             coefH = cxsi2
         else if (nitcfg.eq.1) then
             coefH = 0.
         end if

      else if (b.eq.2) then ! Runge-Kutta ordre 2

         if (nitcfg.eq.2) then
            coefH =-1. 
         else if (nitcfg.eq.1) then 
            coefH = 0.
         end if

      else if (b.eq.1) then ! Runge-Kutta ordre 1

         if (nitcfg.eq.1) then
             coefH = 0.
         end if

      else if (b.eq.5) then       ! Runge-Kutta ordre 3
       
         if (nitcfg.eq.2) then
             coefH = cxsi1_
         else if (nitcfg.eq.3) then
             coefH = cxsi2_
         else if (nitcfg.eq.1) then
             coefH = 0.
         else if (nitcfg.eq.4) then
             coefH = cxsi3_
         else if (nitcfg.eq.5) then
             coefH = cxsi4_
         end if

      else if (b.eq.12) then       ! Runge-Kutta ordre 3
       
         if (nitcfg.eq.1) then
             coefH = 0
         else if (nitcfg.eq.2) then
             coefH = cxsi1_
         else if (nitcfg.eq.3) then
             coefH = cxsi2_
         else if (nitcfg.eq.4) then
             coefH = cxsi3_
         else if (nitcfg.eq.5) then
             coefH = cxsi4_
         else if (nitcfg.eq.6) then
             coefH = cxsi5_
         else if (nitcfg.eq.7) then
             coefH = cxsi6_
         else if (nitcfg.eq.8) then
             coefH = cxsi7_
         else if (nitcfg.eq.9) then
             coefH = cxsi8_
         else if (nitcfg.eq.10) then
             coefH = cxsi9_
         else if (nitcfg.eq.11) then
             coefH = cxsi10_
         else if (nitcfg.eq.12) then
             coefH = cxsi11_
         else if (nitcfg.eq.13) then
             coefH = cxsi12_
         else if (nitcfg.eq.14) then
             coefH = cxsi13_
         end if


      end if

c***********************************************************************
Calcul des coeffcients : Explicit local (RK 2)
c***********************************************************************

      else if (param_int(EXPLOC)==1.and.param_int(RK)==2) then

       if (ind==2) then ! zones de + gd niveau en temps

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         end if

       else  !zones dont le niveau en temps n'est pas le + gd 

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH

         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind/2) then
            coefH=-1.
            !print*, 'coefH=',coefH
         end if   

       end if
      

c***********************************************************************
Calcul des coeffcients : Tang & Warnecke (RK 2)
c***********************************************************************

      else if (param_int(EXPLOC)==5.and.param_int(RK)==2) then

       if (ind==2) then ! zones de + gd niveau en temps

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         end if

       else  !zones dont le niveau en temps n'est pas le + gd 

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH

         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind/2) then
            coefH=-1.
            !print*, 'coefH=',coefH
         end if   

       end if
      

c*************************Constantinescu RK2******************************

      else if (param_int(EXPLOC)==2.and.param_int(RK)==2) then

       if (ind==2) then ! zones de + gd niveau en temps

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.!!!-3 pour les 4 dernieres colonnes !!!
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
         end if


       else  !zones dont le niveau en temps n'est pas le + gd 

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.!!!-3 pour les 4 dernieres colonnes !!!
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind/2) then
            coefH=1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind/2+1) then
            coefH=0.
            !print*, 'coefH=',coefH
         end if   

       end if



c*************************Constantinescu ac rk3***************************

      else if (param_int(EXPLOC)==3.and.param_int(RK)==2) then
         
         if (ind==3) then ! zones de + gd niveau en temps

            if (MOD(nitcfg,ind).eq.0) then
              coefH = cxsi2  
            elseif (MOD(nitcfg,ind).eq.1) then
              coefH = 0.   
            elseif (MOD(nitcfg,ind).eq.2) then
              coefH = cxsi1
            end if
            
       
         else  !zones dont le niveau en temps n'est pas le + gd 

            if (MOD(nitcfg,ind).eq.0) then
              coefH = cxsi2      
            elseif (MOD(nitcfg,ind).eq.1) then
              coefH = 0.                           
            elseif (MOD(nitcfg,ind).eq.2) then
              coefH = cxsi1 
           elseif (MOD(nitcfg,ind).eq.ind/2) then
              coefH = cxsi2  
           elseif (MOD(nitcfg,ind).eq.ind/2+1) then
              coefH = 0.      
           elseif (MOD(nitcfg,ind).eq.ind/2+2) then
              coefH = cxsi1        
           end if

         end if


c*********************** dt local ordre 3*******************************      
      else if (param_int(EXPLOC)==2.and.param_int(RK)==3) then

         if (MOD(nitcfg,ind)==1) then
            coefH=0.
            !print*,'coeff= ',coefH
         elseif(MOD(nitcfg,ind)==ind/2) then
            coefH = cxsi1
            !print*,'coeff= ',coefH
         elseif(MOD(nitcfg,ind)==ind-1) then
            coefH = cxsi2
            !print*,'coeff= ',coefH
         !elseif (MOD(nitcfg,ind)==ind/4.and.ind.ne.4) then
         !   coefH = cxsi3
         !   print*,'coeff= ',coefH
         end if


c*********************RK2 local test************************************


      else if (param_int(EXPLOC)==4.and.param_int(RK)==2) then

       if (ind==2) then ! zones de + gd niveau en temps

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         end if

       else  !zones dont le niveau en temps n'est pas le + gd 

         if (MOD(nitcfg,ind).eq.0) then
            coefH=-1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.1) then
            coefH=0.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind/2) then
            coefH=-1.
            !print*, 'coefH=',coefH
         else if (MOD(nitcfg,ind).eq.ind-1) then
            coefH=0.
            !print*, 'coefH=',coefH
         end if   

       end if
   

      end if

c*****************************************************************************



      IF(param_int(ITYPCP).le.1.or.nitcfg.eq.1.or.b==4) THEN

 
          if(neq.eq.5) then

            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
#include      "FastS/Compute/loopI_begin.for"
                 drodm(l          )=  0.
                 drodm(l +ndimdx  )=  0.
                 drodm(l +ndimdx*2)=  0.
                 drodm(l +ndimdx*3)=  0.
                 drodm(l +ndimdx*4)=  0.
               enddo
            enddo
            enddo

          else

            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
#include      "FastS/Compute/loopI_begin.for"
                 drodm(l          )=  0.
                 drodm(l +ndimdx  )=  0.
                 drodm(l +ndimdx*2)=  0.
                 drodm(l +ndimdx*3)=  0.
                 drodm(l +ndimdx*4)=  0.
                 drodm(l +ndimdx*5)=  0.
               enddo
            enddo
            enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!! Initialisation des cell fictives adjacentes aux frontieres !!!!!!!
        !!!!!!!!!!!!!!!!!(conservativite du schema dt local)!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c$$$        if (param_int(EXPLOC)==2.and.param_int(RK)==3.and.
c$$$     &    param_int(NIJK+4)==0) then
c$$$
c$$$         do ne=1,neq
c$$$           vg = ndimdx*(ne-1)
c$$$               do  j = ind_loop(3), ind_loop(4)
c$$$                 l  = inddm( ind_loop(1)-1 , j, k) 
c$$$                 drodm(l+vg)=0.
c$$$                 l  = inddm( ind_loop(2)+1 , j, k) 
c$$$                 drodm(l+vg)=0.
c$$$                end do
c$$$ 
c$$$
c$$$               do  i = ind_loop(1), ind_loop(2)
c$$$                 l  = inddm( i , ind_loop(3)-1, k) 
c$$$                 drodm(l+vg)=0.
c$$$                 l  = inddm( i , ind_loop(4)+1, k) 
c$$$                 drodm(l+vg)=0.
c$$$                end do
c$$$
c$$$          end do
c$$$
c$$$        else if (param_int(EXPLOC)==2.and.param_int(RK)==3.and.
c$$$     &    param_int(NIJK+4).ne.0) then
c$$$
c$$$         do ne=1,neq
c$$$           vg = ndimdx*(ne-1)
c$$$             do  k = ind_loop(5), ind_loop(6)
c$$$               do  j = ind_loop(3), ind_loop(4)
c$$$                 l  = inddm( ind_loop(1)-1 , j, k) 
c$$$                 drodm(l+vg)=0.
c$$$                 l  = inddm( ind_loop(2)+1 , j, k) 
c$$$                 drodm(l+vg)=0.
c$$$                end do
c$$$              end do
c$$$
c$$$             do  k = ind_loop(5), ind_loop(6)
c$$$               do  i = ind_loop(3), ind_loop(4)
c$$$                 l  = inddm( i , ind_loop(3)-1, k) 
c$$$                 drodm(l+vg)=0.
c$$$                 l  = inddm( 1 , ind_loop(4)+1, k) 
c$$$                 drodm(l+vg)=0.
c$$$                end do
c$$$              end do
c$$$
c$$$             do  j = ind_loop(3), ind_loop(4)
c$$$               do  i = ind_loop(1), ind_loop(2)
c$$$                 l  = inddm( i , j, ind_loop(5)-1) 
c$$$                 drodm(l+vg)=0.
c$$$                 l  = inddm( i , j, ind_loop(6)+1) 
c$$$                 drodm(l+vg)=0.
c$$$                end do
c$$$              end do
c$$$            end do
c$$$          end if



          endif


      ELSE

          if(neq.eq.5) then

            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
#include      "FastS/Compute/loopI_begin.for"

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!Schema de Constantinescu !!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

               if(param_int(EXPLOC)==2.and.param_int(RK)==2.and.
     &              MOD(nitcfg,ind).eq.0.and.ind.ne.b) then
                  
                  if( param_int(LEVEL).lt.param_int(LEVELD)
     &                 .and.param_int(LEVELG).lt.param_int(LEVEL) ) then

                     if(l.ge.lij+1+ind_loop(2)-ind_loop(1)-3) then
                        
                        coefH = -3.
                     else  

                        coefH = -1.
                     end if
                  end if

                  if( param_int(LEVELG).ge.param_int(LEVEL)
     &                 .and.param_int(LEVEL).ge.param_int(LEVELD) ) then

                     if(l.le.lij+4) then          

                        coefH = -3.
                     else 

                        coefH = -1.
                     end if
                  end if

               end if

              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                 drodm(l          )=  drodm(l          )*coefH
                 drodm(l +ndimdx  )=  drodm(l +ndimdx  )*coefH
                 drodm(l +ndimdx*2)=  drodm(l +ndimdx*2)*coefH
                 drodm(l +ndimdx*3)=  drodm(l +ndimdx*3)*coefH
                 drodm(l +ndimdx*4)=  drodm(l +ndimdx*4)*coefH
               enddo
            enddo
            enddo

          else

            do  k = ind_loop(5), ind_loop(6)
            do  j = ind_loop(3), ind_loop(4)
#include      "FastS/Compute/loopI_begin.for"
                 drodm(l          )=  drodm(l          )*coefH
                 drodm(l +ndimdx  )=  drodm(l +ndimdx  )*coefH
                 drodm(l +ndimdx*2)=  drodm(l +ndimdx*2)*coefH
                 drodm(l +ndimdx*3)=  drodm(l +ndimdx*3)*coefH
                 drodm(l +ndimdx*4)=  drodm(l +ndimdx*4)*coefH
                 drodm(l +ndimdx*5)=  drodm(l +ndimdx*5)*coefH
               enddo
            enddo
            enddo

          endif


      ENDIF



      



c$$$      if (param_int(EXPLOC)==2.and.param_int(RK)==3.and.
c$$$     & ind.ne.4 .and. mod(nitcfg,ind)==2) then 
c$$$
c$$$         do  ne= 1, neq
c$$$          vg = ndimdx*(ne-1)
c$$$          do  k = ind_loop(5), ind_loop(6)
c$$$          do  j = ind_loop(3), ind_loop(4)
c$$$
c$$$
c$$$           drodm(l +vg)=  drodmbis(l+vg)
c$$$           drodm(l +vg)=  drodm(l + vg)*coefH
c$$$           
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$        enddo
c$$$
c$$$      end if

        !!!!! Conservatvite du schema dt local d'ordre 3!!!!
c$$$         if (param_int(RK)==3.and.param_int(EXPLOC)==2) then
c$$$
c$$$
c$$$          if (mod(nitcfg,ind)==1) then
c$$$                 do  ne= 1, neq
c$$$                     vg = ndimdx*(ne-1)
c$$$                     do  k = ind_loop(5), ind_loop(6)
c$$$                       do  j = ind_loop(3), ind_loop(4)
c$$$                        l  = inddm( ind_loop(1)-1 , j, k) 
c$$$                        drodm(l+vg)=0.
c$$$                        l  = inddm( ind_loop(2)+1 , j, k) 
c$$$                        drodm(l+vg)=0.
c$$$                        end do
c$$$                      end do
c$$$                    end do
c$$$
c$$$            
c$$$          else if (mod(nitcfg,ind)==ind/2.or.
c$$$     & mod(nitcfg,ind)==ind-1) then
c$$$                 do  ne= 1, neq
c$$$                     vg = ndimdx*(ne-1)
c$$$                     do  k = ind_loop(5), ind_loop(6)
c$$$                       do  j = ind_loop(3), ind_loop(4)
c$$$                        l  = inddm( ind_loop(1)-1 , j, k) 
c$$$                        drodm(l+vg)=coefH*drodm(l+vg)
c$$$                        !drodm(l+vg)=0.
c$$$                        l  = inddm( ind_loop(2)+1 , j, k) 
c$$$                        drodm(l+vg)=coefH*drodm(l+vg)
c$$$                        !drodm(l+vg)=0.
c$$$                        end do
c$$$                      end do
c$$$                    end do
c$$$
c$$$               end if
c$$$          end if




      if (param_int(EXPLOC)==1) then
         ind_loop(2)=ind_loop(2)-1
         ind_loop(1)=ind_loop(1)+1
      end if

      ! if (param_int(EXPLOC)==2 .and. param_int(NIJK+4)==0) then
      !   ind_loop(2)=ind_loop(2)-1
      !   ind_loop(1)=ind_loop(1)+1
      !   ind_loop(4)=ind_loop(4)-1
      !   ind_loop(3)=ind_loop(3)+1
      ! else  if(param_int(EXPLOC)==2 .and. param_int(NIJK+4).ne.0) then
      !   ind_loop(2)=ind_loop(2)-1
      !   ind_loop(1)=ind_loop(1)+1
      !   ind_loop(4)=ind_loop(4)-1
      !   ind_loop(3)=ind_loop(3)+1
      !   ind_loop(6)=ind_loop(6)-1
      !   ind_loop(5)=ind_loop(5)+1
      !end if

      end
