      tn   = 1./max(rot,1.e-27)
      rotx = rotx*tn
      roty = roty*tn
      rotz = rotz*tn

c     calcul des coordonnees des vecteurs de la base cellule
c     renormalisee mais pas reorthogonalisee
c     direction i
       tx=.5*(ti(lt,1)+ti(lt+inci_mtr,1))
       ty=.5*(ti(lt,2)+ti(lt+inci_mtr,2))
       tz=.5*(ti(lt,3)+ti(lt+inci_mtr,3))
c     -- renormalisation
       tn=1./max(sqrt(tx**2+ty**2+tz**2),1.e-27)
       tx=tx*tn
       ty=ty*tn
       tz=tz*tn

       nx=rotx*tx+roty*ty+rotz*tz

c     direction j
       tx=.5*(tj(lt,1)+tj(lt+incj_mtr,1))
       ty=.5*(tj(lt,2)+tj(lt+incj_mtr,2))
       tz=.5*(tj(lt,3)+tj(lt+incj_mtr,3))
c     -- renormalisation
       tn=1./max(sqrt(tx**2+ty**2+tz**2),1.e-27)
       tx=tx*tn
       ty=ty*tn
       tz=tz*tn

       ny=rotx*tx+roty*ty+rotz*tz

c     direction k
       tx=.5*(tk(lt,1)+tk(lt+inck_mtr,1))
       ty=.5*(tk(lt,2)+tk(lt+inck_mtr,2))
       tz=.5*(tk(lt,3)+tk(lt+inck_mtr,3))
c     -- renormalisation
       tn=1./max(sqrt(tx**2+ty**2+tz**2),1.e-27)
       tx=tx*tn
       ty=ty*tn
       tz=tz*tn

       nz=rotx*tx+roty*ty+rotz*tz

c     calcul coordonnees de la direction de rot dans repere cellule
       adelta1 = sqrt(nx*nx*dy*dz+ny*ny*dz*dx+nz*nz*dx*dy)

       adelta1 = (.5+sign(.5,rot-1.e-27))*adelta1
     &          +(.5-sign(.5,rot-1.e-27))*(vol(lvo)**(1./3.))

