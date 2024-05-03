      IF(param_int(ITYPZONE).eq.3) THEN !2D

         if(param_int(NEQ).eq.5) then
            
            do j = ind_loop_lu(3), ind_loop_lu(4)
               do i = ind_loop_lu(1), ind_loop_lu(2)

                  l  =   inddm(i, j, 1)
                  ls = indssor(i, j, 1)

                  ssortmp(ls, 1) = drodm_in(l, 1) + ssor(ls, 1)
                  ssortmp(ls, 2) = drodm_in(l, 2) + ssor(ls, 2)
                  ssortmp(ls, 3) = drodm_in(l, 3) + ssor(ls, 3)
                  ssortmp(ls, 5) = drodm_in(l, 5) + ssor(ls, 5)

                  ssor(ls, 1) = 0.
                  ssor(ls, 2) = 0.
                  ssor(ls, 3) = 0.
                  ssor(ls, 5) = 0.

               enddo
            enddo

         else                   !SA

            do j = ind_loop_lu(3), ind_loop_lu(4)
               do i = ind_loop_lu(1), ind_loop_lu(2)

                  l  =   inddm(i, j, 1)
                  ls = indssor(i, j, 1)

                  ssortmp(ls, 1) = drodm_in(l, 1) + ssor(ls, 1)
                  ssortmp(ls, 2) = drodm_in(l, 2) + ssor(ls, 2)
                  ssortmp(ls, 3) = drodm_in(l, 3) + ssor(ls, 3)
                  ssortmp(ls, 5) = drodm_in(l, 5) + ssor(ls, 5)
                  ssortmp(ls, 6) = drodm_in(l, 6) + ssor(ls, 6)

                  ssor(ls, 1) = 0.
                  ssor(ls, 2) = 0.
                  ssor(ls, 3) = 0.
                  ssor(ls, 5) = 0.
                  ssor(ls, 6) = 0.

               enddo
            enddo

         endif

      else                      !3D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l =   inddm(i, j, k)
                     ls= indssor(i, j, k)

                     ssortmp(ls, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(ls, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(ls, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(ls, 4) = drodm_in(l, 4) + ssor(ls, 4)
                     ssortmp(ls, 5) = drodm_in(l, 5) + ssor(ls, 5)

                     ssor(ls, 1) = 0.
                     ssor(ls, 2) = 0.
                     ssor(ls, 3) = 0.
                     ssor(ls, 4) = 0.
                     ssor(ls, 5) = 0.

                  enddo
               enddo
            enddo

         else                   !SA

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l  =   inddm(i, j, k)
                     ls = indssor(i, j, k)

                     ssortmp(ls, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(ls, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(ls, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(ls, 4) = drodm_in(l, 4) + ssor(ls, 4)
                     ssortmp(ls, 5) = drodm_in(l, 5) + ssor(ls, 5)
                     ssortmp(ls, 6) = drodm_in(l, 6) + ssor(ls, 6)

                     ssor(ls, 1) = 0.
                     ssor(ls, 2) = 0.
                     ssor(ls, 3) = 0.
                     ssor(ls, 4) = 0.
                     ssor(ls, 5) = 0.
                     ssor(ls, 6) = 0.

                  enddo
               enddo
            enddo

         endif
      endif
