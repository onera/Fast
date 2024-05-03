      amut   =rop(l,1)*anutild*fvv1
      xmuprov=max(amut,0.)
      xmut(l)=  amulam
     &        + min(xmuprov,param_real(SA_REAL+SA_RATIOM-1)*amulam)
