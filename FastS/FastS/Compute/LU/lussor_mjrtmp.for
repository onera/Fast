      IF(param_int(ITYPZONE).eq.3) THEN !2D

         if(param_int(NEQ).eq.5) then
            
            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l = inddm(i, j, k)
                     ls = indssor(i, j, k, i_size, j_size)

                     ssortmp(l, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(l, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(l, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(l, 5) = drodm_in(l, 5) + ssor(ls, 5)

                     ssor(ls, 1) = 0.
                     ssor(ls, 2) = 0.
                     ssor(ls, 3) = 0.
                     ssor(ls, 5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l = inddm(i, j, k)
                     ls = indssor(i, j, k, i_size, j_size)

                     ssortmp(l, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(l, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(l, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(l, 5) = drodm_in(l, 5) + ssor(ls, 5)
                     ssortmp(l, 6) = drodm_in(l, 6) + ssor(ls, 6)

                     ssor(ls, 1) = 0.
                     ssor(ls, 2) = 0.
                     ssor(ls, 3) = 0.
                     ssor(ls, 5) = 0.
                     ssor(ls, 6) = 0.

                  enddo
               enddo
            enddo

         endif

      else !3D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l = inddm(i, j, k)
                     ls = indssor(i, j, k, i_size, j_size)

                     ssortmp(l, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(l, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(l, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(l, 4) = drodm_in(l, 4) + ssor(ls, 4)
                     ssortmp(l, 5) = drodm_in(l, 5) + ssor(ls, 5)

                     ssor(ls, 1) = 0.
                     ssor(ls, 2) = 0.
                     ssor(ls, 3) = 0.
                     ssor(ls, 4) = 0.
                     ssor(ls, 5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     l = inddm(i, j, k)
                     ls = indssor(i, j, k, i_size, j_size)

                     ssortmp(l, 1) = drodm_in(l, 1) + ssor(ls, 1)
                     ssortmp(l, 2) = drodm_in(l, 2) + ssor(ls, 2)
                     ssortmp(l, 3) = drodm_in(l, 3) + ssor(ls, 3)
                     ssortmp(l, 4) = drodm_in(l, 4) + ssor(ls, 4)
                     ssortmp(l, 5) = drodm_in(l, 5) + ssor(ls, 5)
                     ssortmp(l, 6) = drodm_in(l, 6) + ssor(ls, 6)

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
