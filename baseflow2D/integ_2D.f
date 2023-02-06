      program integ 20140617

      ! this program reads the output of the baseflow (based_*.bin) and
      ! calculates the values of displacement thickness (delta1),
      ! momentum thickness (delta2) and shape factor H12 = delta1 / delta2
      implicit none
      include '../par.for'
      character*15 nome
      integer i, j, my_rank, shift
      real*8 delta(imax), delta1(imax), delta2(imax), H12(imax),
     &       ux(ptsx,jmax), uy(ptsx,jmax), wz(ptsx,jmax),
     &       uxt(imax,jmax), uyt(imax,jmax), wzt(imax,jmax),
     &       var1(imax,jmax), var2(imax,jmax),
     &       a, b, c, det, fc1, fc2, umax, x, y(jmax)
      
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'based_',my_rank,'.bin'
        open(2, file = nome, form = 'unformatted')
        read(2) ux, uy, wz
        close (unit = 2)

        shift = my_rank * (ptsx - inter - 1)
        do j = 1, jmax
          do i = 1, ptsx
            uxt(i+shift,j) = ux(i,j)
            uyt(i+shift,j) = uy(i,j)
            wzt(i+shift,j) = wz(i,j)
          end do
        end do
      end do

      do j = 1, jmax
        if (stf .eq. 1.d0) then 
          ! without stretching
          y(j) = dy * dble(j-1)
         else
          ! with stretching
          y(j) = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if        
      end do

      do i = 1, imax
        umax = uxt(i,jmax)
        delta(i) = 0.d0
        do j = 1, jmax
          var1(i,j) = 1.d0 - uxt(i,j) / umax
          var2(i,j) = uxt(i,j) / umax * (1.d0 - uxt(i,j) / umax)
          if (uxt(i,j)/umax .lt. 0.99d0) then
            delta(i) = y(j)
          end if
        end do
      end do

      do i = 1, imax
        fc1 = 0.d0
        fc2 = 0.d0
        do j = 2, jmax - 1, 2
          det = y(j-1)**2 * ( y(j) - y(j+1) )
     &        + y(j-1) * ( y(j+1)**2 - y(j)**2 )
     &        + y(j) * y(j+1) * ( y(j) - y(j+1) )
          a   = (1.d0 / det) * ( var1(i,j-1) * (y(j)   - y(j+1))
     &                         + var1(i,j)   * (y(j+1) - y(j-1))
     &                         + var1(i,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var1(i,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var1(i,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var1(i,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det)
     &        * ( var1(i,j-1) * y(j)   * y(j+1) * (y(j)   - y(j+1))
     &          + var1(i,j)   * y(j-1) * y(j+1) * (y(j+1) - y(j-1))
     &          + var1(i,j+1) * y(j)   * y(j-1) * (y(j-1) - y(j)))
          fc1 = fc1 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &              + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &              +  c         * (y(j+1)    - y(j-1))

          a   = (1.d0 / det) * ( var2(i,j-1) * (y(j)   - y(j+1))
     &                         + var2(i,j)   * (y(j+1) - y(j-1))
     &                         + var2(i,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var2(i,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var2(i,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var2(i,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det)
     &        * ( var2(i,j-1) * y(j)   * y(j+1) * (y(j)   - y(j+1))
     &          + var2(i,j)   * y(j-1) * y(j+1) * (y(j+1) - y(j-1))
     &          + var2(i,j+1) * y(j)   * y(j-1) * (y(j-1) - y(j)))
          fc2 = fc2 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3)
     &              + (b / 2.d0) * (y(j+1)**2 - y(j-1)**2)
     &              +  c         * (y(j+1)    - y(j-1))
        end do
        if (my_form .eq. 0) then
          delta1(i) = fc1 * L_1 * 1.d3 / dsqrt(fac_y) ! 1000 -> mm
          delta2(i) = fc2 * L_1 * 1.d3 / dsqrt(fac_y) ! 1000 -> mm
         else
          delta1(i) = fc1 / dsqrt(fac_y)
          delta2(i) = fc2 / dsqrt(fac_y)
        end if
        H12(i) = delta1(i) / delta2(i)
      end do

      open (1, file = 'H12.dat',status = 'unknown')
      write(1,*) 'VARIABLES="x","delta","delta1","delta2","H12"'
      write(1,*) 'ZONE I=',imax,',  F=POINT'
      write(*,*) ' The results are stored in the file H12.dat'
      do i = 1, imax
        x = x0 + dble(i-1) * dx
        write(1,5) x, delta(i), delta1(i), delta2(i), H12(i)
      end do
      close (unit = 1)
 
      if (my_form .eq. 0) then
        open (2, file = 'delta1.dat',status = 'unknown')
        do i = 1, imax
          write(2,6) delta1(i)
        end do
        close (unit = 2)
      end if
 
    5 format(1x, 5d25.17)
    6 format(1x, 1d25.17)

      return
      end
