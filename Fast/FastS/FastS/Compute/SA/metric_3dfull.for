               !calcul de delta et gradient vitesse
               tix  = ti(lt,1)
               tiy  = ti(lt,2)
               tiz  = ti(lt,3)
               tix1 = ti(lt+inci_mtr,1)
               tiy1 = ti(lt+inci_mtr,2)
               tiz1 = ti(lt+inci_mtr,3)

               tci  =       sqrt( tix*tix+tiy*tiy+tiz*tiz)
               tci  = tci + sqrt( tix1*tix1+tiy1*tiy1+tiz1*tiz1) 
               tci  = max(.5*tci , 1.e-27)
               dx   = vol(lvo)/tci
c
               tjx  = tj(lt,1)
               tjy  = tj(lt,2)
               tjz  = tj(lt,3)
               tjx1 = tj(lt+incj_mtr,1)
               tjy1 = tj(lt+incj_mtr,2)
               tjz1 = tj(lt+incj_mtr,3)
c
               tcj  =       sqrt(  tjx*tjx + tjy*tjy + tjz*tjz)
               tcj  = tcj + sqrt( tjx1*tjx1+tjy1*tjy1+tjz1*tjz1) 
               tcj  = max(.5*tcj,1.e-27)
               dy   = vol(lvo)/tcj

               tkx  = tk(lt,1)
               tky  = tk(lt,2)
               tkz  = tk(lt,3)
               tkx1 = tk(lt+inck_mtr,1)
               tky1 = tk(lt+inck_mtr,2)
               tkz1 = tk(lt+inck_mtr,3)


               tck  =       sqrt(  tkx*tkx + tky*tky + tkz*tkz)
               tck  = tck + sqrt( tkx1*tkx1+tky1*tky1+tkz1*tkz1) 
               tck  = max(.5*tck,1.e-27)
               dz   = vol(lvo)/tck
               
               sph2 = max(dx,dy,dz)
 
