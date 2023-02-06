cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                    subroutines transforms binary data in asc                 c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 20140617

      implicit none
      include '../par.for'
      character*20 nm
      integer i, j, my_rank, shift
      real*8 x(imax), y(imax,jmax), xx, yy, delta1(imax),
     &       uxt(imax,jmax), uyt(imax,jmax), 
     &       wzt(imax,jmax), tht(imax,jmax),
     &        ux(ptsx,jmax),  uy(ptsx,jmax),  
     &        wz(ptsx,jmax),  th(ptsx,jmax),
     &       S_z(ptsx,jmax), S_zt(ptsx,jmax)

      ! 0 -> Zero Curvature
      ! 1 -> Gortler vortices simulations without heat transfer
      ! 2 -> Gortler vortices simulations with heat transfers
      select case (my_form)
        
        case(1) ! -> Gortler vortices simulations withour heat transfer
          ! Baseflow variables data
          do my_rank = 0, np - 1
            write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
            open(1, file = nm, form = 'unformatted')
            read(1) ux, uy, wz
            close (unit = 1)
            shift = my_rank * (ptsx - inter - 1)
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j) = ux(i,j)
                uyt(i+shift,j) = uy(i,j)
                wzt(i+shift,j) = wz(i,j)
              end do
            end do
          end do

          ! writes data to spatial space to be open by tecplot
          open (3, file = 'basens.dat',status = 'unknown')
          write(3,*) 'VARIABLES="x","y","ux","uy","wz"'
          write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretching
              yy = dble(j-1) * dy
             else
              ! with stretching
              yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1, imax
              xx = x0 + dble(i-1) * dx
              write(3,5) xx, yy, uxt(i,j), uyt(i,j), wzt(i,j)
            end do
          end do
          close (unit = 3)

          open(1, file = 'basens.bin', form = 'unformatted')
          write(1) uxt, uyt, wzt
          close (unit = 1)


        case(2) ! -> Gortler vortices simulations with heat transfer
          ! Baseflow variables data
          do my_rank = 0, np - 1
            write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
            open(1, file = nm, form = 'unformatted')
            read(1) ux, uy, wz, th
            close (unit = 1)
            shift = my_rank * (ptsx - inter - 1)
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j) = ux(i,j)
                uyt(i+shift,j) = uy(i,j)
                wzt(i+shift,j) = wz(i,j)
                tht(i+shift,j) = th(i,j)
              end do
            end do
          end do

          ! writes data to spatial space to be open by tecplot
          open (3, file = 'basens.dat',status = 'unknown')
          write(3,*) 'VARIABLES="x","y","ux","uy","wz","theta"'
          write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

          do j = 1, jmax
            if(stf .eq. 1.d0) then
              ! without stretching
              yy = dble(j-1) * dy
             else
              ! with stretching
              yy = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1, imax
              xx = x0 + dble(i-1) * dx
              write(3,6) xx, yy, uxt(i,j), uyt(i,j), wzt(i,j), tht(i,j)
            end do
          end do
          close (unit = 3)

          open(1, file = 'basens.bin', form = 'unformatted')
          write(1) uxt, uyt, wzt, tht
          close (unit = 1)

        end select

    5 format(1x, 2d14.6, 3e21.9e3)
    6 format(1x, 2d14.6, 4e21.9e3)
        
      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                 end subroutines transforms binary data in asc                c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
