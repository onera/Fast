      tn = 1./max(rot,1.e-27)
      nx = rotx*tn
      ny = roty*tn
      nz = rotz*tn

c     calcul coordonnees de la direction de rot dans repere cellule
      adelta1 = sqrt(nx*nx*dy*dz+ny*ny*dz*dx+nz*nz*dx*dy)

       adelta1 = (.5+sign(.5,rot-1.e-27))*adelta1
     &          +(.5-sign(.5,rot-1.e-27))*(vol(lvo)**(1./3.))


