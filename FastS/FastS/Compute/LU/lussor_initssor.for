      IF(param_int(ITYPZONE).eq.3) THEN !2D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5), ind_loop_lu(6)
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v5) = 0.
                     ssor(l + v6) = 0.

                  enddo
               enddo
            enddo

         endif

      else !3D

         if(param_int(NEQ).eq.5) then

            do k = ind_loop_lu(5) - 1, ind_loop_lu(6) + 1
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v4) = 0.
                     ssor(l + v5) = 0.

                  enddo
               enddo
            enddo

         else !SA

            do k = ind_loop_lu(5) - 1, ind_loop_lu(6) + 1
               do j = ind_loop_lu(3) - 1, ind_loop_lu(4) + 1

                  lij = indssor(ind_loop_lu(1), j, k)

                  do l = lij-1, lij + ind_loop_lu(2) - ind_loop_lu(1)+1

                     ssor(l) = 0.
                     ssor(l + v2) = 0.
                     ssor(l + v3) = 0.
                     ssor(l + v4) = 0.
                     ssor(l + v5) = 0.
                     ssor(l + v6) = 0.

                  enddo
               enddo
            enddo

         endif
      endif
