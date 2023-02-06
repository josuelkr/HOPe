cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                  subroutines transforms binary data in asc                   c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program formated 04092012

      implicit none
      include '../par.for'
      character*20 nm
      integer i, j, my_rank, shift
      real*8 x, y, ux(imax,jmax), uy(imax,jmax), wz(imax,jmax),
     &       th(imax,jmax)

      select case (my_form)
        case(1) ! -> Gortler vortices simulation without heat transfer

          open(1, file = 'base_fs.bin', form = 'unformatted')
          read(1) ux, uy, wz
          close (unit = 1)
          ! writes data to spacial space to be open by tecplot
          open (3, file = 'base_fs.dat',status = 'unknown')
          write(3,*) 'VARIABLES="x","y","velu","vely","vortz"'
          write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'
        
          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretchting
              y = dble(j-1) * dy
            else
              ! with stretchting
              y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1, imax
              x = x0 + dble(i-1) * dx
              write(3,5) x, y, ux(i,j), uy(i,j), wz(i,j)
            end do
          end do
        
        case(2) ! -> Gortler vortices simulation with heat transfer

          open(1, file = 'base_fs.bin', form = 'unformatted')
          read(1) ux, uy, wz, th
          close (unit = 1)
        
          ! writes data to spacial space to be open by tecplot
          open (3, file = 'base_fs.dat',status = 'unknown')
          write(3,*) 'VARIABLES="x","y","velu","vely","vortz","theta"'
          write(3,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretchting
              y = dble(j-1) * dy
            else
              ! with stretchting
              y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do i = 1, imax
              x = x0 + dble(i-1) * dx
              write(3,6) x, y, ux(i,j), uy(i,j), wz(i,j), th(i,j)
            end do
          end do
      end select
 
    5 format(1x, 2d14.6, 3e21.9e3)
    6 format(1x, 2d14.6, 4e21.9e3)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                end subroutines transforms binary data in asc                 c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
