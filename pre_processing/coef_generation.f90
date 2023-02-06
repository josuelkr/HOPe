!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

program coef_generation

 use solvers
 use derivative

 implicit none

 call init_base_flow

end program coef_generation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

subroutine init_base_flow

  use derivative
  use atwall_calculation
  use solvers

  implicit none
  integer      :: i, j, k, m_ptsx
  real(kind=8) :: rhs(imax), dudx(imax,jmax), dvdy(imax,jmax), x, y, dudx2(imax,jmax)
  real(kind=8) :: dvdy2(imax,jmax), sp(7,lvls), cp(8,lvls), pp(4,lvls), lp(5,lvls)
  character(len=17) :: nm

  call at_wall_coef_generation
  call derivative_coefficient_generation
  do i = 1 , lvls
   call poisson_coefficient_generation(stf**2**(i-1))
   sp(:,i) = sp_poi_coef
   cp(:,i) = cp_poi_coef
   pp(:,i) = pp_poi_coef
   lp(:,i) = lp_poi_coef
  end do
  open(1, file = 'coefs.bin', form = 'unformatted')
   write(1) fp_fd_coef_e
   write(1) sp_fd_coef
   write(1) cp_fd_coef
   write(1) pp_fd_coef
   write(1) lp_fd_coef
   write(1) fp_sd_coef
   write(1) sp_sd_coef
   write(1) cp_sd_coef
   write(1) pp_sd_coef
   write(1) lp_sd_coef
   write(1) sp
   write(1) cp
   write(1) pp
   write(1) lp
   write(1) w_at_w_coef
   write(1) dwydy_coef
   write(1) sp_integ_coef
   write(1) mp_integ_coef
   write(1) pp_integ_coef
   write(1) lp_integ_coef
  close(unit = 1)
 
  write(*,*) "Done!"

end subroutine init_base_flow

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

