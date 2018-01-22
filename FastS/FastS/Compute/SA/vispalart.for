c***********************************************************************
c     $Date: 2013-02-04 19:27:20 +0100 (lun. 04 févr. 2013) $
c     $Revision: 38 $
c     $Author: IvanMary $
c***********************************************************************
      subroutine vispalart(ndom,  param_int, param_real, ind_loop,
     &                     ti,tj,tk,vol,dlng, xmut,rop)
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

      INTEGER_E ndom,ind_loop(6), param_int(0:*)

      REAL_E xmut( param_int(NDIMDX) )
      REAL_E  rop( param_int(NDIMDX) , param_int(NEQ) )

      REAL_E ti( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tj( param_int(NDIMDX_MTR) , param_int(NEQ_IJ) ),
     &       tk( param_int(NDIMDX_MTR) , param_int(NEQ_K ) )
      REAL_E  vol( param_int(NDIMDX_MTR) )
      REAL_E dlng( param_int(NDIMDX) )

      REAL_E param_real(0:*)
c Var loc 
      INTEGER_E incmax,i,j,k,l,ltij,lij,lt,lvo

      REAL_E amulam,anulam,fv1,fvv1,anutild,amut,xmuprov,
     &  chi,ad1, s,temp01,cmus1,coesut,t1_1

#include "FastS/formule_mtr_param.h"
#include "FastS/formule_param.h"

c.....formulation originelle
      fv1(s)     = (s**3)/(s**3+SA_CV1)
 
      cmus1  =    param_real(VISCO+4)
      temp01 = 1./param_real(VISCO+3)
      coesut =    param_real(VISCO+2) * (1.+cmus1*temp01)

      !
      !
      !Modele Spalart
      !
      !
      IF (param_int(SA_INT + SA_IDES-1 ).le.1) THEN  !SA

#include     "FastS/Compute/loop_begin.for" 
       
#include     "FastS/Compute/mulam.for" 
#include     "FastS/Compute/SA/chi.for"
              fvv1   = fv1(chi)
#include     "FastS/Compute/SA/xmut.for" 

#include     "FastS/Compute/loop_end.for" 

      ENDIF


      end
      !
      !
      !Modeles ZDES1,2,3,....
      !
      !
C      ELSE
C 
C        !!idist = des_int(SA_DIST)
C        !! idist= 1: c1d =1, c2d =0, c3d = 0
C        !! idist= 2: c1d =0, c2d =1, c3d = 0
C        !! idist= 3: c1d =0, c2d =0, c3d = 1
C
C        if( param_int(SA_INT+ SA_IZGRIS-1 ).eq.0) then
C
C         call viszdes_izgris0(ndom, param_int, param_real, ind_loop,
C     &                        ti,tj,tk,vol,dlng, xmut,rop)
C
C        else
C
C          !Si izgris=1, alors idist n'a pas d'influence: delta1 = dlng
C#include  "FastS/Compute/loop_begin.for" 
C
C#include     "FastS/Compute/mulam.for" 
C#include     "FastS/Compute/SA/chi.for" 
C             ad1     = max(dlng(l), 1.e-27)  
C#include     "SA/fvv1_izgris1.for" 
C#include     "FastS/Compute/SA/xmut.for" 
C
C#include  "FastS/Compute/loop_end.for"  
C
C        endif

