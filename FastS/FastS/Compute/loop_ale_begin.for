#ifndef E_SCALAR_COMPUTER
CDIR$ IVDEP
!CDIR NODEP
CVD$  NODEPCHK
             do l=incmax+1,ndimdx-incmax
#else
             do k = ind_loop(5), ind_loop(6)
             do j = ind_loop(3), ind_loop(4)
                lij  =       inddm( ind_loop(1) , j, k)-1
                ltij = lij - indmtr(ind_loop(1) , j, k)+1
                lvij = lij - indven(ind_loop(1) , j, k)+1
!DEC$ IVDEP
                do l = lij+1, lij+1 +  ind_loop(2)- ind_loop(1)
#endif 
                    lt = l  - ltij
                    lv = l  - lvij
                    lvo= lt
