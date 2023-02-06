      program integ 20121123

      ! this program reads the output of fs_bound.f and calculates the
      ! values of displacement thickness (alfa1), momentum thickness (alfa2)
      ! and calculates the H12=alfa1/alfa2
      implicit none
      include '../par.for'
      integer i, j
      real*8 delta1(imax), delta2(imax), H12(imax), fc1, fc2, umax,
     &       ux(imax,jmax), uy(imax,jmax), wz(imax,jmax), y(jmax),
     &       a, b, c, det, var1(imax,jmax), var2(imax,jmax), x
      
      do j = 1, jmax
        if (stf .ne. 1.d0) then 
         y(j) = dy * (stf**(j-1)-1.d0) / (stf-1.d0)
        else
         y(j) = dy * dble(j-1)
        endif        
      enddo

      open(1,file='basens.bin',form='unformatted')
      read(1) ux, uy, wz
      close(unit=1)
      
      do i = 1, imax
        umax = ux(i,jmax)
        do j = 1, jmax
          var1(i,j) = 1.d0 - ux(i,j) / umax
          var2(i,j) = ux(i,j) / umax * (1.d0 - ux(i,j) / umax)
        end do
      end do

      do i = 1, imax
        fc1 = 0.d0
        fc2 = 0.d0
        do j = 2, jmax - 1, 2
          det = y(j-1)**2 * ( y(j) - y(j+1) )                      
     &        + y(j-1) * ( y(j+1)**2 - y(j)**2)               
     &        + y(j) * y(j+1) * ( y(j)-y(j+1) ) 
          a   = (1.d0 / det) * ( var1(i,j-1) * (y(j)   - y(j+1))
     &                         + var1(i,j)   * (y(j+1) - y(j-1))
     &                         + var1(i,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var1(i,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var1(i,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var1(i,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det) * ( var1(i,j-1) *
     &                           y(j)  * y(j+1) * (y(j)   - y(j+1))
     &                         + var1(i,j)   *
     &                           y(j-1)* y(j+1) * (y(j+1) - y(j-1))
     &                         + var1(i,j+1) *
     &                           y(j)  * y(j-1) * (y(j-1) - y(j)))
          fc1 = fc1 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3) + 
     &                (b / 2.d0) * (y(j+1)**2 - y(j-1)**2) + 
     &                 c         * (y(j+1)    - y(j-1))

          a   = (1.d0 / det) * ( var2(i,j-1) * (y(j)   - y(j+1))
     &                         + var2(i,j)   * (y(j+1) - y(j-1))
     &                         + var2(i,j+1) * (y(j-1) - y(j)))
          b   = (1.d0 / det) * ( var2(i,j-1) * (y(j+1)**2 - y(j)**2)
     &                         + var2(i,j)   * (y(j-1)**2 - y(j+1)**2)
     &                         - var2(i,j+1) * (y(j-1)**2 - y(j)**2))
          c   = (1.d0 / det) * ( var2(i,j-1) *
     &                           y(j)  * y(j+1) * (y(j)   - y(j+1))
     &                         + var2(i,j)   *
     &                           y(j-1)* y(j+1) * (y(j+1) - y(j-1))
     &                         + var2(i,j+1) *
     &                           y(j)  * y(j-1) * (y(j-1) - y(j)))
          fc2 = fc2 + (a / 3.d0) * (y(j+1)**3 - y(j-1)**3) + 
     &                (b / 2.d0) * (y(j+1)**2 - y(j-1)**2) + 
     &                 c         * (y(j+1)    - y(j-1))
        end do
        delta1(i) = fc1 / dsqrt(fac_y)
        delta2(i) = fc2 / dsqrt(fac_y)
        H12(i)    = delta1(i) / delta2(i)
      end do

      open (1, file = 'H12.dat',status = 'unknown')
      write(1,*) 'VARIABLES="x","delta1","delta2","H12"'
      write(1,*) 'ZONE I=',imax,',  F=POINT'
      write(*,*) ' The results are stored in the file H12.dat' 
      do i = 1, imax
        x = (dble(i-1) * dx + x0)
        write(1,5) x, delta1(i), delta2(i), H12(i)
      end do
      close (unit=1)
 
    5 format(1x,4d25.17)
      return
      end
