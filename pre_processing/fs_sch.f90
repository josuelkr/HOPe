
 module fs

  use constants

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 contains

  subroutine fpp_finder(fpp, eta_adp, i)

   implicit none
   real(kind=8), intent(inout) :: fpp, eta_adp
   integer, intent(in)         :: i
   integer                     :: j, cont, steps_eta
   real(kind=8)                :: f_eta(6), delta_fpp, ls_err

   ls_err = 1.d0
   cont   = 1

   do while ( (dabs(ls_err) .gt. gtol_fs) .and. (cont .lt. 10000) )
     delta_fpp = 100.d0
     steps_eta = eta_adp / deta + 1

     do while (dabs(delta_fpp/fpp) .gt. tol_fs)
       f_eta(1) = 0.d0
       f_eta(2) = 0.d0
       f_eta(3) = fpp
       f_eta(4) = 0.d0
       f_eta(5) = 0.d0
       f_eta(6) = 1.d0
      
       do j = 2, steps_eta
         call rk4(f_eta, deta, i)
       end do

       delta_fpp = ( - f_eta(5) * f_eta(2) - ( f_eta(6) * f_eta(3) ) + f_eta(5) ) / ( f_eta(5)**2 + f_eta(6)**2 )
       fpp       = fpp + delta_fpp
     end do

     ls_err  = ( 1.d0 - f_eta(2) )**2 + f_eta(3)**2
     eta_adp = eta_adp + eta_variation
     cont    = cont + 1
   end do

   eta_adp = eta_adp - eta_variation
   print *, cont-1, beta_fs(i), fpp, eta_adp, ls_err

  end subroutine fpp_finder

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine baseflow_fs(fpp, thp)

   implicit none
   real(kind=8), intent(in) :: fpp(imax), thp(imax)
   real(kind=8)             :: f_eta(6), h, m, ue, xad, yad, eta, deta0, x, y
   real(kind=8)             :: ux(imax,jmax), uy(imax,jmax), wz(imax,jmax), th(imax,jmax)
   integer                  :: i, j, pt, k

   do i = 1, imax
     m        = beta_fs(i) / ( 2.d0 - beta_fs(i) )
     xad      = x0 + dble(i-1) * dx
     ue       = xad ** m
     eta      = 0.d0
     deta0    = ( dy0 / dsqrt(fac_y) ) * dsqrt(0.5d0 * Re * ue * (m + 1.d0) / xad)

     f_eta(1) = 0.d0
     f_eta(2) = 0.d0
     f_eta(3) = fpp(i)
     f_eta(4) = 0.d0
     f_eta(5) = thp(i)
     f_eta(6) = 0.d0

     ux(i,1)  = 0.d0
     uy(i,1)  = 0.d0
     wz(i,1)  = ue * deta0 * f_eta(3) / ( dy0 / dsqrt(fac_y) )
     wz(i,1)  = wz(i,1) / dsqrt(fac_y)   ! adimensionalization
     th(i,1)  = 0.d0

     do j = 2, jmax
       h  = deta0 * (stf**(j-2))
       pt = int( h  / 1.d-4 )
       h  = h / dble(pt)

       do k = 1, pt
         call rk4_theta(f_eta, h, i)
       end do

       if (stf .ne. 1.d0) then
         eta = deta0 * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        else
         eta = deta0 * (j-1)
       end if

       yad     = eta / dsqrt(0.5d0 * Re * ue * (m + 1.d0) / xad)
       ux(i,j) = ue * f_eta(2)
       uy(i,j) = - dsqrt( 0.5d0 * ue* ( m + 1.d0 ) / (xad * Re)) * f_eta(1) &
                 - 0.5d0 * ( m - 1.d0 ) * yad * ue * f_eta(2) / xad
       uy(i,j) = uy(i,j) * dsqrt(fac_y)
       wz(i,j) = ue * dsqrt(0.5d0 * (m + 1.d0) * Re * ue / xad) * f_eta(3) !                              &
!              + 0.125d0 * (m - 1.d0)**2 * yad**2 * dsqrt(2.d0 * (m + 1.d0) / (Re * xad)) * f_eta(3) / xad**2  &
!              + 0.25d0 * (3.d0 * m**2 - 4.d0 * m + 1.d0) * yad / (xad**2 * dsqrt(Re)) * f_eta(2)              &
!              + 0.25d0 * (m - 1.d0) * dsqrt(2.d0 * (m + 1.d0) /(Re * xad**3)) * f_eta(1)
       wz(i,j) = wz(i,j) / dsqrt(fac_y)   ! adimensionalization
       th(i,j) = f_eta(4)
     end do
   end do

   if (my_form .eq. 2) then 
     open(1, file = 'base_fs.bin', form = 'unformatted')
     write(1) ux, uy, wz, th
     close(unit = 1)
     
!    open (2, file = 'base_fs.dat', status = 'unknown')
!    write(2,*) 'VARIABLES="x","y","ux","uy","wz","theta"'
!    write(2,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'
!    
!    do j = 1, jmax
!      if (stf .eq. 1.d0) then
!        ! without stretching
!        y = dble(j-1) * dy0
!       else
!        ! with stretching
!        y = dy0 * (stf**(j-1) - 1.d0) / (stf - 1.d0)
!      end if
!      do i = 1, imax
!        x = x0 + dble(i-1) * dx
!        write(2,'(1x, 2d14.6, 4d24.16)') x, y, ux(i,j), uy(i,j), wz(i,j), th(i,j)
!      end do
!    end do
!    close(unit = 2)
    else
     open(1, file = 'base_fs.bin', form = 'unformatted')
     write(1) ux, uy, wz
     close(unit = 1)

!    open (2, file = 'base_fs.dat',status = 'unknown')
!    write(2,*) 'VARIABLES="x","y","ux","uy","wz"'
!    write(2,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'

!    do j = 1, jmax
!      if (stf .eq. 1.d0) then
!        ! without stretching
!        y = dble(j-1) * dy0
!       else
!        ! with stretching
!        y = dy0 * (stf**(j-1) - 1.d0) / (stf - 1.d0)
!      end if
!      do i = 1, imax
!        x = x0 + dble(i-1) * dx
!        write(2,*) x, y, ux(i,j), uy(i,j), wz(i,j)
!      end do
!    end do
!    close(unit = 2)
   end if
 
  end subroutine baseflow_fs
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
  subroutine delta_calculation(fpp, eta_end, i)

   implicit none
   real(kind=8), intent(in) :: fpp, eta_end
   integer, intent(in)      :: i
   real(kind=8)             :: f_eta(6), h, delta
   real(kind=8)             :: ux, ux_new, variation, delta_eta
   integer                  :: j, Neta

   delta_eta = 1.d-4
   Neta      = int( eta_end  / delta_eta) + 1

   f_eta(1) = 0.d0
   f_eta(2) = 0.d0
   f_eta(3) = fpp
   f_eta(4) = 0.d0
   f_eta(5) = 0.d0
   f_eta(6) = 1.d0

   ux = 0.d0

   do j = 2, Neta
     call rk4(f_eta, delta_eta, i)
     if (f_eta(2) .ge. 0.99d0) then
       delta = dble(j-1) * delta_eta
       write(*,*) '99% = ', delta
       return
     end if
   end do
 
  end subroutine delta_calculation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine rk4(f_eta, h, i)

   implicit none
   integer, intent(in)         :: i
   real(kind=8), intent(inout) :: f_eta(6)
   real(kind=8), intent(in)    :: h
   real(kind=8)                :: hh, h6, dfd_eta(6), f_etat(6), dfd_eta2(6)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs(f_eta, dfd_eta, i)
   f_etat  = f_eta   + hh   * dfd_eta

   call derivs(f_etat, dfd_eta2, i) 
   f_etat  = f_eta   + hh   * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs(f_etat, dfd_eta2, i)
   f_etat  = f_eta   + h    * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs(f_etat, dfd_eta2, i)
   f_eta   = f_eta   + h6   * ( dfd_eta + dfd_eta2 )

  end subroutine rk4

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine derivs (f_eta, dfd_eta, i)

   implicit none
   real(kind=8), intent(in)  :: f_eta(6)
   real(kind=8), intent(out) :: dfd_eta(6)
   integer, intent(in)       :: i

   ! Eq. (9.8) Schlichting page 164
   dfd_eta(1) =   f_eta(2)
   dfd_eta(2) =   f_eta(3)
   dfd_eta(3) = - f_eta(1) * f_eta(3) - beta_fs(i) * (1.d0 - f_eta(2)**2)
   dfd_eta(4) =   f_eta(5)
   dfd_eta(5) =   f_eta(6)
   dfd_eta(6) = - f_eta(4) * f_eta(3) - f_eta(1) * f_eta(6) + 2.d0 * beta_fs(i) * f_eta(2) * f_eta(5)

  end subroutine derivs

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine rk42(f_eta, h, i)

   implicit none
   integer, intent(in)         :: i
   real(kind=8), intent(inout) :: f_eta(6)
   real(kind=8), intent(in)    :: h
   real(kind=8)                :: hh, h6, dfd_eta(6), f_etat(6), dfd_eta2(6)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs2(f_eta, dfd_eta, i)
   f_etat   = f_eta + hh * dfd_eta

   call derivs2(f_etat, dfd_eta2, i) 
   f_etat   = f_eta + hh * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs2(f_etat, dfd_eta2, i)
   f_etat   = f_eta + h * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs2(f_etat, dfd_eta2, i)
   f_eta   = f_eta + h6 * ( dfd_eta + dfd_eta2 )
 
  end subroutine rk42
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine derivs2 (f_eta, dfd_eta, i)

   implicit none
   real(kind=8), intent(in)  :: f_eta(6)
   real(kind=8), intent(out) :: dfd_eta(6)
   integer, intent(in)       :: i
   real(kind=8)              :: m

   m = beta_fs(i) / (2.d0 - beta_fs(i))   ! pag 164 do livro do Schlichting, m significa o parametro de adimensionalização pressão/gradiente.

   dfd_eta(1) =   f_eta(2)
   dfd_eta(2) =   f_eta(3)
   dfd_eta(3) = - f_eta(1) * f_eta(3) * 0.5d0 * ( m + 1.d0 ) - m * ( 1.d0 - f_eta(2) * f_eta(2) )
   dfd_eta(4) =   f_eta(5)
   dfd_eta(5) =   f_eta(6)
   dfd_eta(6) =   ( - f_eta(4) * f_eta(3) - f_eta(1) * f_eta(6)) + (2.d0 * beta_fs(i) * f_eta(2) * f_eta(5))

  end subroutine derivs2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! THETA DISTRIBUTION CALCULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine thetap_finder(thp, fpp, i)

   implicit none
   integer, intent(in)       :: i
   real(kind=8), intent(in)  :: fpp
   real(kind=8), intent(out) :: thp
   integer                   :: j
   real(kind=8)              :: x1, x2, xacc, f, fh, fl, temp, xh, xl, df, dxold, dxx

   x1   = 0.0d0
   x2   = 1.5d0
!  x1   = fpp / Pr - 0.25d0
!  x2   = fpp / Pr + 0.25d0
   xacc = 5.d-16

   call funcd_th(x1, fl, df, fpp, i)
   call funcd_th(x2, fh, df, fpp, i)

   if (fl .lt. 0.d0) then
     xl = x1
     xh = x2
    else
     xh = x1
     xl = x2
   end if
   thp   = 0.5d0 * (x1 + x2)
   dxold = dabs(x2 - x1)
   dxx   = dxold
   call funcd_th(thp, f, df, fpp, i)
   do j = 1, 50
     if ( (((thp-xh)*df-f)*((thp-xl)*df-f) .ge. 0.d0) .or. &
          (dabs(2.d0*f) .gt. dabs(dxold*df)) ) then
       dxold = dxx
       dxx   = 0.5d0 * (xh - xl)
       thp   = xl + dxx
      else
       dxold = dxx
       dxx   = f / df
       temp  = thp
       thp   = thp - dxx
     end if

     call funcd_th(thp, f, df, fpp, i)
     if (dabs(f) .lt. xacc) return
     if (f .lt. 0.d0) then
       xl = thp
      else
       xh = thp
     end if

   end do

   write(*,*)'thetap_finder exceeding maximum iteractions -> thp'

  end subroutine thetap_finder

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine funcd_th(x, fc, dfc, fpp, i)

   implicit none
   integer, intent(in)       :: i
   real(kind=8), intent(in)  :: x, fpp
   real(kind=8), intent(out) :: fc, dfc
   integer                   :: j, npt
   real(kind=8)              :: f_eta(6), dxl, fc2, h

   h    = 1.d-4
   npt  = 150000

   f_eta(1) = 0.d0
   f_eta(2) = 0.d0
   f_eta(3) = fpp
   f_eta(4) = 0.d0
   f_eta(5) = x
   f_eta(6) = 0.d0
   do j = 2, npt
     call rk4_theta(f_eta, h, i)
   end do
   fc = f_eta(4) - 1.d0
 
   dxl = 1.d-6

   f_eta(1) = 0.d0
   f_eta(2) = 0.d0
   f_eta(3) = fpp
   f_eta(4) = 0.d0
   f_eta(5) = x + dxl
   f_eta(6) = 0.d0
   do j = 2, npt
     call rk4_theta(f_eta, h, i)
   end do
   fc2 = f_eta(4) - 1.d0
 
   dfc = (fc2 - fc) / dxl
 
  end subroutine funcd_th

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine rk4_theta(f_eta, h, i)

   implicit none
   integer, intent(in)         :: i
   real(kind=8), intent(inout) :: f_eta(6)
   real(kind=8), intent(in)    :: h
   real(kind=8)                :: hh, h6, dfd_eta(6), f_etat(6), dfd_eta2(6)

   hh = h * 0.5d0
   h6 = h / 6.d0

   call derivs_theta(f_eta, dfd_eta, i)
   f_etat  = f_eta + hh * dfd_eta

   call derivs_theta(f_etat, dfd_eta2, i)
   f_etat  = f_eta + hh * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs_theta(f_etat, dfd_eta2, i)
   f_etat  = f_eta + h * dfd_eta2
   dfd_eta = dfd_eta + 2.d0 * dfd_eta2

   call derivs_theta(f_etat, dfd_eta2, i)
   f_eta   = f_eta + h6 * ( dfd_eta + dfd_eta2 )

  end subroutine rk4_theta

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

  subroutine derivs_theta(f_eta, dfd_eta, i)

   implicit none
   real(kind=8), intent(in)  :: f_eta(6)
   real(kind=8), intent(out) :: dfd_eta(6)
   integer, intent(in)       :: i

   ! Eq. (9.8) Schlichting page 164
   dfd_eta(1) =   f_eta(2)
   dfd_eta(2) =   f_eta(3)
   dfd_eta(3) = - f_eta(1) * f_eta(3) - (beta_fs(i) * (1.d0 - f_eta(2)**2))
   dfd_eta(4) =   f_eta(5)
   dfd_eta(5) = - f_eta(1) * f_eta(5) * Pr
   dfd_eta(6) =   0.d0

  end subroutine derivs_theta

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! END OF THETA DISTRIBUTION CALCULATION
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 end module fs
