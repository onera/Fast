               !calcul de delta et gradient vitesse
               tix  = ti(lt,1)
               tiy  = ti(lt,2)
               tix1 = ti(lt+inci_mtr,1)
               tiy1 = ti(lt+inci_mtr,2)

               tci  =       sqrt( tix*tix+tiy*tiy)
               tci  = tci + sqrt( tix1*tix1+tiy1*tiy1) 
               tci  = max(.5*tci , 1.e-27)
               dx   = vol(lvo)/tci
c
               tjx  = tj(lt,1)
               tjy  = tj(lt,2)
               tjx1 = tj(lt+incj_mtr,1)
               tjy1 = tj(lt+incj_mtr,2)
c
               tcj  =       sqrt(  tjx*tjx + tjy*tjy )
               tcj  = tcj + sqrt( tjx1*tjx1+tjy1*tjy1) 
               tcj  = max(.5*tcj,1.e-27)
               dy   = vol(lvo)/tcj

               sph2 = max(dx,dy)
 
