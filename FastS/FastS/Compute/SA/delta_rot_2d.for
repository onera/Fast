      tn   = 1./max(rot,1.e-27)
      nz   = rotz*tn

c     calcul coordonnees de la direction de rot dans repere cellule
       adelta1 = sqrt(nz*nz*dx*dy)

       adelta1 = (.5+sign(.5,rot-1.e-27))*adelta1
     &          +(.5-sign(.5,rot-1.e-27))*(vol(lvo)**(1./3.))

