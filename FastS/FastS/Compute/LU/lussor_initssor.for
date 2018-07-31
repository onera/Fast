      IF(param_int(ITYPZONE).eq.3) THEN !2D

         if(param_int(NEQ).eq.5) then

            do j = ind_loop_lu(3), ind_loop_lu(4)
               do i = ind_loop_lu(1), ind_loop_lu(2)

                  ls = indssor(i, j, 1, i_size, j_size)

                  ssor(ls, 1) = 0.
                  ssor(ls, 2) = 0.
                  ssor(ls, 3) = 0.
                  ssor(ls, 5) = 0.

               enddo
            enddo

         else !SA

            do j = ind_loop_lu(3), ind_loop_lu(4)
               do i = ind_loop_lu(1), ind_loop_lu(2)

                  ls = indssor(i, j, 1, i_size, j_size)

                  ssor(ls, 1) = 0.
                  ssor(ls, 2) = 0.
                  ssor(ls, 3) = 0.
                  ssor(ls, 5) = 0.
                  ssor(ls, 6) = 0.

               enddo
            enddo

         endif

      else !3D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3), ind_loop_lu(4)
                  do i = ind_loop_lu(1), ind_loop_lu(2)

                     ls = indssor(i, j, k, i_size, j_size)

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

                     ls = indssor(i, j, k, i_size, j_size)

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
