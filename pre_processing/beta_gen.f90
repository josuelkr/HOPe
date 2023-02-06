 program beta_gen

  use constants
  real(kind = 8) :: x
  integer        :: i

  ! writes beta_fs in a file
  open (1, file = '../beta_fs.dist',status = 'unknown')
  do i = 1, imax
   x = (i - 1.d0) * dx
!  beta_fs(i) = 3.9156626510d-4 * dble(i) -  6.03915662651d-2 
!  beta_fs(i) = -0.06d0
   beta_fs(i) = 0.0d0
  write(1,'(1d17.9)') beta_fs(i)
  end do
  close(unit = 1)
  
 end program beta_gen
