c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine viszdes_izgris0(ndom, param_int, param_real, ind_loop, 
     &                           ti,tj,tk,vol,dlng, xmut,rop)  
c***********************************************************************
c_P                          O N E R A
c     ACT
c_A    Calcul de la contribution des termes sources 
c_A    pour l'equation de Spalart Allmaras 
c
c     VAL
c_V    Gaz parfait mono-espece
c_V    Navier-Stokes
c
c     INP
c_I    ndom      : numero du domaine calcule
c_I    dvardc(6) : gradients de nutild primitif aux centres
c_I                des cellules
c_I    dvardc(5) : gradients de ronutild conservatif aux centres
c_I                des cellules
c_I    vol       : volumes
c_I    rotn      : norme du rotationnel
c_I    dlng      : distance a la paroi
c
c     OUT
c
c     I/O
c_I    drodm    : terme source de l'equation de Spalart Allmaras 
c***********************************************************************
      implicit  none

#include "FastS/param_solver.h"

      INTEGER_E ndom,ind_loop(6),param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E dlng( param_int(NDIMDX) )
      REAL_E  rop( param_int(NDIMDX) , param_int(NEQ) )

      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )
      REAL_E vol(param_int(NDIMDX_MTR))

      REAL_E param_real(0:*)
c Var loc 
      INTEGER_E incmax,i,j,k,l,icoe_pos,inci,incj,inck,
     &  l1,l2,l3,l4,l5,l6,ltij,lij,lt,inci_mtr,
     & incj_mtr,inck_mtr, lvo

      REAL_E amulam,anulam,fv1,fvv1,ad1,adelta1,adelta2,adelta3,hwn,
     &  amutild,anutild,amut,xmuprov,tci,tcj,tck,
     &  chi,r,r2,g,fwg1,fwg,c1,cw1,
     & f1,f2,s,t,
     & temp01,cmus1,coesut,t1,t1_1,voldes,
     & tjx, tjy,tjz,tjx1,tjy1,tjz1,si,sj,sk,
     & tix, tiy,tiz,tix1,tiy1,tiz1,sph2,
     & tkx, tky,tkz,tkx1,tky1,tkz1

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

c.....formulation originelle
      fv1(s)     = (s**3)/(s**3+SA_CV1)

      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

      c1     = SA_CB2/SA_SIGMA
      cw1    = (SA_CB1/SA_CKARM/SA_CKARM)+(1.+SA_CB2)/SA_SIGMA
 
      inci_mtr = param_int(NIJK_MTR)
      incj_mtr = param_int(NIJK_MTR+1)
      inck_mtr = param_int(NIJK_MTR+2)


      IF(param_int(SA_INT+ SA_IDIST-1).eq.1) then  

#include  "FastS/Compute/loop_begin.for" 

#include     "FastS/Compute/mulam.for" 
#include     "FastS/Compute/SA/chi.for" 
             ad1     = max(dlng(l), 1.e-27)  
#include     "FastS/Compute/SA/delta1.for" 
#include     "FastS/Compute/SA/fvv1_izgris0.for" 
#include     "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for" 


      ELSEIF(param_int(SA_INT+ SA_IDIST-1).eq.2) then  

        if( param_int(ITYPZONE).eq.0) then

#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_3d.for" 
              adelta1= adelta2 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for" 

        elseif( param_int(ITYPZONE).eq.1) then

#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_3dhomo.for" 
              adelta1= adelta2 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for"


        elseif( param_int(ITYPZONE).eq.2) then

          !metric
          lt  = indmtr(1 , 1, 1)
          lvo = lt
          tix = ti(lt,1)
          tjy = tj(lt,1)
          tkz = tk(lt,1)
          si      = abs (tix)
          sj      = abs (tjy)
          sk      = abs (tkz)
          sph2    = min(si,sj,sk)
          sph2    = vol(lvo)/sph2
          voldes  = vol(lvo) 
          adelta1 = (voldes)**(1./3.)
          adelta2 = max(sph2,1.e-27)

          !Verif
#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
              adelta1= adelta2 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for"

        else
#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_2d.for" 
              adelta1= adelta2 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for" 
        endif

      ELSE   !param_int(SA_INT+ SA_IDIST-1)=3

        if( param_int(ITYPZONE).eq.0) then

#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_3d.for" 
              !attention provisoire cas particulier PP avec normale en y
#include      "FastS/Compute/SA/delta3.for" 
              adelta1= adelta3 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for" 

        elseif( param_int(ITYPZONE).eq.1) then

#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_3dhomo.for" 
              !attention provisoire cas particulier PP avec normale en y
#include      "FastS/Compute/SA/delta3.for" 
              adelta1= adelta3 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for"


        elseif( param_int(ITYPZONE).eq.2) then

          !metric
          lt  = indmtr(1 , 1, 1)
          lvo = lt
          tix = ti(lt,1)
          tjy = tj(lt,1)
          tkz = tk(lt,1)
          si      = abs (tix)
          sj      = abs (tjy)
          sk      = abs (tkz)
          sph2    = min(si,sj,sk)
          sph2    = vol(lvo)/sph2
          voldes  = vol(lvo) 
          adelta1 = (voldes)**(1./3.)
          adelta2 = max(sph2,1.e-27)

          !attention provisoire cas particulier PP avec normale en y
          tcj=sj

#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
              !attention provisoire cas particulier PP avec normale en y
#include      "FastS/Compute/SA/delta3.for" 
              adelta1= adelta3 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for"

        else
#include  "FastS/Compute/loop_begin.for" 

#include      "FastS/Compute/mulam.for" 
#include      "FastS/Compute/SA/chi.for" 
              ad1     = max(dlng(l), 1.e-27)  
#include      "FastS/Compute/SA/delta1.for" 
#include      "FastS/Compute/SA/delta2_2d.for" 
              !attention provisoire cas particulier PP avec normale en y
#include      "FastS/Compute/SA/delta3.for" 
              adelta1= adelta3 
#include      "FastS/Compute/SA/fvv1_izgris0.for" 
#include      "FastS/Compute/SA/xmut.for" 

#include  "FastS/Compute/loop_end.for" 
        endif

      ENDIF

      end
