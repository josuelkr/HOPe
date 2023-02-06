
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

module derivative

 use derivative_commom_variables
 use solvers

 contains

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine first_point_first_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hp3, hp4, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 = dy
  hp2 = hp1 + stf    * dy
  hp3 = hp2 + stf**2 * dy
  hp4 = hp3 + stf**3 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0

  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) = - hp1
  m(2,5) = - hp2
  m(2,6) = - hp3
  m(2,7) = - hp4

  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
  m(3,5) = - 1.d0
  m(3,6) = - 1.d0
  m(3,7) = - 1.d0
 
  m(4,1) =   0.d0
  m(4,2) =   hp1
  m(4,3) =   0.d0
  m(4,4) = - hp1**2 / 2.d0
  m(4,5) = - hp2**2 / 2.d0
  m(4,6) = - hp3**2 / 2.d0
  m(4,7) = - hp4**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) =   hp1**2 / 2.d0
  m(5,3) =   0.d0
  m(5,4) = - hp1**3 / 6.d0
  m(5,5) = - hp2**3 / 6.d0
  m(5,6) = - hp3**3 / 6.d0
  m(5,7) = - hp4**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) =   hp1**3 / 6.d0
  m(6,3) =   0.d0
  m(6,4) = - hp1**4 / 24.d0
  m(6,5) = - hp2**4 / 24.d0
  m(6,6) = - hp3**4 / 24.d0
  m(6,7) = - hp4**4 / 24.d0

  m(7,1) =   0.d0
  m(7,2) =   hp1**4 / 24.d0
  m(7,3) =   0.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) = - hp2**5 / 120.d0
  m(7,6) = - hp3**5 / 120.d0
  m(7,7) = - hp4**5 / 120.d0

 return

 end subroutine first_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_first_point_first_derivative_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_first_point_first_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine second_point_first_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm1, hp1, hp2, hp3, hp4, dy
 
  dy = 1.d0
  ! 'normalized' length
  hm1 = dy
  hp1 =       stf    * dy
  hp2 = hp1 + stf**2 * dy
  hp3 = hp2 + stf**3 * dy
  hp4 = hp3 + stf**4 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  m(1,8) =   0.d0
  m(1,9) =   0.d0

  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   1.d0
  m(2,4) =   hm1
  m(2,5) =   0.d0
  m(2,6) = - hp1
  m(2,7) = - hp2
  m(2,8) = - hp3
  m(2,9) = - hp4
 
  m(3,1) = - hm1
  m(3,2) =   0.d0
  m(3,3) =   hp1
  m(3,4) = - hm1**2 / 2.d0
  m(3,5) =   0.d0
  m(3,6) = - hp1**2 / 2.d0
  m(3,7) = - hp2**2 / 2.d0
  m(3,8) = - hp3**2 / 2.d0
  m(3,9) = - hp4**2 / 2.d0
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) = - 1.d0
  m(4,5) = - 1.d0
  m(4,6) = - 1.d0
  m(4,7) = - 1.d0
  m(4,8) = - 1.d0
  m(4,9) = - 1.d0
 
  m(5,1) =   hm1**2 / 2.d0
  m(5,2) =   0.d0
  m(5,3) =   hp1**2 / 2.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   0.d0
  m(5,6) = - hp1**3 / 6.d0
  m(5,7) = - hp2**3 / 6.d0
  m(5,8) = - hp3**3 / 6.d0
  m(5,9) = - hp4**3 / 6.d0
 
  m(6,1) = - hm1**3 / 6.d0
  m(6,2) =   0.d0
  m(6,3) =   hp1**3 / 6.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hp1**4 / 24.d0
  m(6,7) = - hp2**4 / 24.d0
  m(6,8) = - hp3**4 / 24.d0
  m(6,9) = - hp4**4 / 24.d0

  m(7,1) =   hm1**4 / 24.d0
  m(7,2) =   0.d0
  m(7,3) =   hp1**4 / 24.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) = - hp1**5 / 120.d0
  m(7,7) = - hp2**5 / 120.d0
  m(7,8) = - hp3**5 / 120.d0
  m(7,9) = - hp4**5 / 120.d0

  m(8,1) = - hm1**5 / 120.d0
  m(8,2) =   0.d0
  m(8,3) =   hp1**5 / 120.d0
  m(8,4) = - hm1**6 / 720.d0
  m(8,5) =   0.d0
  m(8,6) = - hp1**6 / 720.d0
  m(8,7) = - hp2**6 / 720.d0
  m(8,8) = - hp3**6 / 720.d0
  m(8,9) = - hp4**6 / 720.d0

! m(9,1) =   hm1**6 / 720.d0
! m(9,2) =   0.d0
! m(9,3) =   hp1**6 / 720.d0
! m(9,4) = - hm1**7 / 5040.d0
! m(9,5) =   0.d0
! m(9,6) =   hp1**7 / 5040.d0
! m(9,7) =   hp2**7 / 5040.d0
! m(9,8) =   hp3**7 / 5040.d0
! m(9,9) =   hp4**7 / 5040.d0

  m(9,1) =   0.d0
  m(9,2) =   0.d0
  m(9,3) =   0.d0
  m(9,4) =   0.d0
  m(9,5) = - 1.d0
  m(9,6) =   0.d0
  m(9,7) =   0.d0
  m(9,8) =   0.d0
  m(9,9) =   0.d0
 
  return

 end subroutine second_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_second_point_first_derivative_matrix_coef(b, n)

  implicit none
  ! right-hand side of first derivative coef.
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
! do i = 2, n
  do i = 1, n
   b(i) = 0.d0
  end do
  b(1) = 1.d0
! b(2) = 1.d0
  b(9) = 2.5d0
 
  return

 end subroutine rhs_second_point_first_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine center_point_first_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8)              :: hp1, hp2, hm1, hm2
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: dy

  dy = 1.d0  
  ! 'normalized' length
  hm1 =       stf
  hm2 = dy  + stf
  hp1 =       stf**2
  hp2 = hp1 + stf**3

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   1.d0
  m(1,3) =   1.d0
  m(1,4) =   hm2
  m(1,5) =   hm1
  m(1,6) =   0.d0
  m(1,7) = - hp1
  m(1,8) = - hp2
 
  m(2,1) =   0.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) =   0.d0
  m(2,5) =   0.d0
  m(2,6) =   0.d0
  m(2,7) =   0.d0
  m(2,8) =   0.d0

  m(3,1) = - hm1
  m(3,2) =   0.d0
  m(3,3) =   hp1
  m(3,4) = - hm2**2 / 2.d0
  m(3,5) = - hm1**2 / 2.d0
  m(3,6) =   0.d0
  m(3,7) = - hp1**2 / 2.d0
  m(3,8) = - hp2**2 / 2.d0

  m(4,1) =   hm1**2 / 2.d0
  m(4,2) =   0.d0
  m(4,3) =   hp1**2 / 2.d0
  m(4,4) =   hm2**3 / 6.d0
  m(4,5) =   hm1**3 / 6.d0
  m(4,6) =   0.d0
  m(4,7) = - hp1**3 / 6.d0
  m(4,8) = - hp2**3 / 6.d0
 
  m(5,1) = - hm1**3 / 6.d0
  m(5,2) =   0.d0
  m(5,3) =   hp1**3 / 6.d0
  m(5,4) = - hm2**4 / 24.d0
  m(5,5) = - hm1**4 / 24.d0
  m(5,6) =   0.d0
  m(5,7) = - hp1**4 / 24.d0
  m(5,8) = - hp2**4 / 24.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - 1.d0
  m(6,5) = - 1.d0
  m(6,6) = - 1.d0
  m(6,7) = - 1.d0
  m(6,8) = - 1.d0

  m(7,1) =   hm1**4 / 24.d0
  m(7,2) =   0.d0
  m(7,3) =   hp1**4 / 24.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) = - hp1**5 / 120.d0
  m(7,8) = - hp2**5 / 120.d0

  m(8,1) = - hm1**5 / 120.d0
  m(8,2) =   0.d0
  m(8,3) =   hp1**5 / 120.d0
  m(8,4) = - hm2**6 / 720.d0
  m(8,5) = - hm1**6 / 720.d0
  m(8,6) =   0.d0
  m(8,7) = - hp1**6 / 720.d0
  m(8,8) = - hp2**6 / 720.d0
 
  return

 end subroutine center_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_center_point_first_derivative_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 3, n
    b(i) = 0.d0
  end do
  b(1) = 0.d0
  b(2) = 1.d0
 
  return

 end subroutine rhs_center_point_first_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine penultimate_point_first_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hm1, hm2, hm3, hm4, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 =       stf**4
  hm1 =       stf**3
  hm2 = hm1 + stf**2
  hm3 = hm2 + stf
  hm4 = hm3 + dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  m(1,8) =   0.d0
  m(1,9) =   0.d0

  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   1.d0
  m(2,4) = - hp1
  m(2,5) =   0.d0
  m(2,6) =   hm1
  m(2,7) =   hm2
  m(2,8) =   hm3
  m(2,9) =   hm4

  m(3,1) =   hp1
  m(3,2) =   0.d0
  m(3,3) = - hm1
  m(3,4) = - hp1**2 / 2.d0
  m(3,5) =   0.d0
  m(3,6) = - hm1**2 / 2.d0
  m(3,7) = - hm2**2 / 2.d0
  m(3,8) = - hm3**2 / 2.d0
  m(3,9) = - hm4**2 / 2.d0
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) = - 1.d0
  m(4,5) = - 1.d0
  m(4,6) = - 1.d0
  m(4,7) = - 1.d0
  m(4,8) = - 1.d0
  m(4,9) = - 1.d0
 
  m(5,1) =   hp1**2 / 2.d0
  m(5,2) =   0.d0
  m(5,3) =   hm1**2 / 2.d0
  m(5,4) = - hp1**3 / 6.d0
  m(5,5) =   0.d0
  m(5,6) =   hm1**3 / 6.d0
  m(5,7) =   hm2**3 / 6.d0
  m(5,8) =   hm3**3 / 6.d0
  m(5,9) =   hm4**3 / 6.d0
 
  m(6,1) =   hp1**3 / 6.d0
  m(6,2) =   0.d0
  m(6,3) = - hm1**3 / 6.d0
  m(6,4) = - hp1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hm1**4 / 24.d0
  m(6,7) = - hm2**4 / 24.d0
  m(6,8) = - hm3**4 / 24.d0
  m(6,9) = - hm4**4 / 24.d0

  m(7,1) =   hp1**4 / 24.d0
  m(7,2) =   0.d0
  m(7,3) =   hm1**4 / 24.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) =   hm1**5 / 120.d0
  m(7,7) =   hm2**5 / 120.d0
  m(7,8) =   hm3**5 / 120.d0
  m(7,9) =   hm4**5 / 120.d0

  m(8,1) =   hp1**5 / 120.d0
  m(8,2) =   0.d0
  m(8,3) = - hm1**5 / 120.d0
  m(8,4) = - hp1**6 / 720.d0
  m(8,5) =   0.d0
  m(8,6) = - hm1**6 / 720.d0
  m(8,7) = - hm2**6 / 720.d0
  m(8,8) = - hm3**6 / 720.d0
  m(8,9) = - hm4**6 / 720.d0

! m(9,1) =   hp1**6 / 720.d0
! m(9,2) =   0.d0
! m(9,3) =   hm1**6 / 720.d0
! m(9,4) =   hp1**7 / 5040.d0
! m(9,5) =   0.d0
! m(9,6) = - hm1**7 / 5040.d0
! m(9,7) = - hm2**7 / 5040.d0
! m(9,8) = - hm3**7 / 5040.d0
! m(9,9) = - hm4**7 / 5040.d0

  m(9,1) =   0.d0
  m(9,2) =   0.d0
  m(9,3) =   0.d0
  m(9,4) =   0.d0
  m(9,5) =   1.d0
  m(9,6) =   0.d0
  m(9,7) =   0.d0
  m(9,8) =   0.d0
  m(9,9) =   0.d0

  return

 end subroutine penultimate_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_penultimate_point_first_derivative_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 2, n
   b(i) = 0.d0
  end do
  b(1) = 1.d0
  b(9) = 2.5d0
 
  return

 end subroutine rhs_penultimate_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine last_point_first_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm1, hm2, hm3, hm4, dy

  dy  = 1.d0
  ! 'normalized' length
  hm1 =       stf**3
  hm2 = hm1 + stf**2
  hm3 = hm2 + stf
  hm4 = hm3 + dy
 
 ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0

  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) =   hm1
  m(2,5) =   hm2
  m(2,6) =   hm3
  m(2,7) =   hm4
  
  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
  m(3,5) = - 1.d0
  m(3,6) = - 1.d0
  m(3,7) = - 1.d0

  m(4,1) =   0.d0
  m(4,2) = - hm1
  m(4,3) =   0.d0
  m(4,4) = - hm1**2 / 2.d0
  m(4,5) = - hm2**2 / 2.d0
  m(4,6) = - hm3**2 / 2.d0
  m(4,7) = - hm4**2 / 2.d0

  m(5,1) =   0.d0
  m(5,2) =   hm1**2 / 2.d0
  m(5,3) =   0.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   hm2**3 / 6.d0
  m(5,6) =   hm3**3 / 6.d0
  m(5,7) =   hm4**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hm1**3 / 6.d0
  m(6,3) =   0.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) = - hm2**4 / 24.d0
  m(6,6) = - hm3**4 / 24.d0
  m(6,7) = - hm4**4 / 24.d0

  m(7,1) =   0.d0
  m(7,2) =   hm1**4 / 24.d0
  m(7,3) =   0.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   hm2**5 / 120.d0
  m(7,6) =   hm3**5 / 120.d0
  m(7,7) =   hm4**5 / 120.d0

  return

 end subroutine last_point_first_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_last_point_first_derivative_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)

  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0

  return

  end subroutine rhs_last_point_first_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine first_point_second_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hp3, hp4, hp5, dy

  dy = 1.d0
  ! 'normalized' length
  hp1 = dy
  hp2 = hp1 + stf    * dy
  hp3 = hp2 + stf**2 * dy
  hp4 = hp3 + stf**3 * dy
  hp5 = hp4 + stf**4 * dy

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  m(1,8) =   0.d0
 
  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) = - hp1**2 / 2.d0
  m(2,5) = - hp2**2 / 2.d0
  m(2,6) = - hp3**2 / 2.d0
  m(2,7) = - hp4**2 / 2.d0
  m(2,8) = - hp5**2 / 2.d0
 
  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
  m(3,5) = - 1.d0
  m(3,6) = - 1.d0
  m(3,7) = - 1.d0
  m(3,8) = - 1.d0

  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - hp1
  m(6,5) = - hp2
  m(6,6) = - hp3
  m(6,7) = - hp4
  m(6,8) = - hp5
 
  m(4,1) =   0.d0
  m(4,2) =   hp1
  m(4,3) =   0.d0
  m(4,4) = - hp1**3 / 6.d0
  m(4,5) = - hp2**3 / 6.d0
  m(4,6) = - hp3**3 / 6.d0
  m(4,7) = - hp4**3 / 6.d0
  m(4,8) = - hp5**3 / 6.d0

  m(5,1) =   0.d0
  m(5,2) =   hp1**2 / 2.d0
  m(5,3) =   0.d0
  m(5,4) = - hp1**4 / 24.d0
  m(5,5) = - hp2**4 / 24.d0
  m(5,6) = - hp3**4 / 24.d0
  m(5,7) = - hp4**4 / 24.d0
  m(5,8) = - hp5**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hp1**3 / 6.d0
  m(7,3) =   0.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) = - hp2**5 / 120.d0
  m(7,6) = - hp3**5 / 120.d0
  m(7,7) = - hp4**5 / 120.d0
  m(7,8) = - hp5**5 / 120.d0

  m(8,1) =   0.d0
  m(8,2) =   hp1**4 / 24.d0
  m(8,3) =   0.d0
  m(8,4) = - hp1**6 / 720.d0
  m(8,5) = - hp2**6 / 720.d0
  m(8,6) = - hp3**6 / 720.d0
  m(8,7) = - hp4**6 / 720.d0
  m(8,8) = - hp5**6 / 720.d0

 return

 end subroutine first_point_second_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_first_point_second_derivative_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_first_point_second_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine second_point_second_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm1, hp1, hp2, hp3, hp4, hp5, dy
 
  dy = 1.d0
  ! 'normalized' length
  hm1 = dy
  hp1 =       stf    * dy
  hp2 = hp1 + stf**2 * dy
  hp3 = hp2 + stf**3 * dy
  hp4 = hp3 + stf**4 * dy
  hp5 = hp4 + stf**5 * dy
 
  ! matrix definition
  m(1,1)  =   1.d0
  m(1,2)  =   0.d0
  m(1,3)  =   0.d0
  m(1,4)  =   0.d0
  m(1,5)  =   0.d0
  m(1,6)  =   0.d0
  m(1,7)  =   0.d0
  m(1,8)  =   0.d0
  m(1,9)  =   0.d0
  m(1,10) =   0.d0
 
  m(2,1)  =   1.d0
  m(2,2)  =   1.d0
  m(2,3)  =   1.d0
  m(2,4)  = - hm1**2 / 2.d0
  m(2,5)  =   0.d0
  m(2,6)  = - hp1**2 / 2.d0
  m(2,7)  = - hp2**2 / 2.d0
  m(2,8)  = - hp3**2 / 2.d0
  m(2,9)  = - hp4**2 / 2.d0
  m(2,10) = - hp5**2 / 2.d0
 
  m(3,1)  = - hm1
  m(3,2)  =   0.d0
  m(3,3)  =   hp1
  m(3,4)  =   hm1**3 / 6.d0
  m(3,5)  =   0.d0
  m(3,6)  = - hp1**3 / 6.d0
  m(3,7)  = - hp2**3 / 6.d0
  m(3,8)  = - hp3**3 / 6.d0
  m(3,9)  = - hp4**3 / 6.d0
  m(3,10) = - hp5**3 / 6.d0
 
  m(4,1)  =   hm1**2 / 2.d0
  m(4,2)  =   0.d0
  m(4,3)  =   hp1**2 / 2.d0
  m(4,4)  = - hm1**4 / 24.d0
  m(4,5)  =   0.d0
  m(4,6)  = - hp1**4 / 24.d0
  m(4,7)  = - hp2**4 / 24.d0
  m(4,8)  = - hp3**4 / 24.d0
  m(4,9)  = - hp4**4 / 24.d0
  m(4,10) = - hp5**4 / 24.d0
 
  m(5,1)  =   0.d0
  m(5,2)  =   0.d0
  m(5,3)  =   0.d0
  m(5,4)  = - 1.d0
  m(5,5)  = - 1.d0
  m(5,6)  = - 1.d0
  m(5,7)  = - 1.d0
  m(5,8)  = - 1.d0
  m(5,9)  = - 1.d0
  m(5,10) = - 1.d0

  m(6,1)  =   0.d0
  m(6,2)  =   0.d0
  m(6,3)  =   0.d0
  m(6,4)  =   0.d0
  m(6,5)  =   0.d0
  m(6,6)  = - 1.d0
  m(6,7)  =   0.d0
  m(6,8)  =   0.d0
  m(6,9)  =   0.d0
  m(6,10) =   0.d0
 
  m(7,1)  = - hm1**3 / 6.d0
  m(7,2)  =   0.d0
  m(7,3)  =   hp1**3 / 6.d0
  m(7,4)  =   hm1**5 / 120.d0
  m(7,5)  =   0.d0
  m(7,6)  = - hp1**5 / 120.d0
  m(7,7)  = - hp2**5 / 120.d0
  m(7,8)  = - hp3**5 / 120.d0
  m(7,9)  = - hp4**5 / 120.d0
  m(7,10) = - hp5**5 / 120.d0
 
  m(8,1)  =   hm1**4 / 24.d0
  m(8,2)  =   0.d0
  m(8,3)  =   hp1**4 / 24.d0
  m(8,4)  = - hm1**6 / 720.d0
  m(8,5)  =   0.d0
  m(8,6)  = - hp1**6 / 720.d0
  m(8,7)  = - hp2**6 / 720.d0
  m(8,8)  = - hp3**6 / 720.d0
  m(8,9)  = - hp4**6 / 720.d0
  m(8,10) = - hp5**6 / 720.d0
 
  m(9,1)  = - hm1**5 / 120.d0
  m(9,2)  =   0.d0
  m(9,3)  =   hp1**5 / 120.d0
  m(9,4)  =   hm1**7 / 5040.d0
  m(9,5)  =   0.d0
  m(9,6)  = - hp1**7 / 5040.d0
  m(9,7)  = - hp2**7 / 5040.d0
  m(9,8)  = - hp3**7 / 5040.d0
  m(9,9)  = - hp4**7 / 5040.d0
  m(9,10) = - hp5**7 / 5040.d0
 
  m(10,1)  =   0.d0
  m(10,2)  =   0.d0
  m(10,3)  =   0.d0
  m(10,4)  =   hm1
  m(10,5)  =   0.d0
  m(10,6)  = - hp1
  m(10,7)  = - hp2
  m(10,8)  = - hp3
  m(10,9)  = - hp4
  m(10,10) = - hp5
 
  return

 end subroutine second_point_second_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_second_point_second_derivative_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
   b(i) = 0.d0
  end do
  b(1) =   1.d0
  b(6) = - 5.25d0
 
  return

 end subroutine rhs_second_point_second_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine center_point_second_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8)              :: hp1, hp2, hm1, hm2
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: dy

  dy = 1.d0  
  ! 'normalized' length
  hm1 =       stf
  hm2 = dy  + stf
  hp1 =       stf**2
  hp2 = hp1 + stf**3

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   1.d0
  m(1,3) =   1.d0
  m(1,4) = - hm2**2 / 2.d0
  m(1,5) = - hm1**2 / 2.d0
  m(1,6) = - 0.d0
  m(1,7) = - hp1**2 / 2.d0
  m(1,8) = - hp2**2 / 2.d0
 
  m(2,1) =   0.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) =   0.d0
  m(2,5) =   0.d0
  m(2,6) =   0.d0
  m(2,7) =   0.d0
  m(2,8) =   0.d0
 
  m(3,1) = - hm1
  m(3,2) =   0.d0
  m(3,3) =   hp1
  m(3,4) =   hm2**3 / 6.d0
  m(3,5) =   hm1**3 / 6.d0
  m(3,6) =   0.d0
  m(3,7) = - hp1**3 / 6.d0
  m(3,8) = - hp2**3 / 6.d0
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) =   hm2
  m(4,5) =   hm1
  m(4,6) =   0.d0
  m(4,7) = - hp1
  m(4,8) = - hp2
 
  m(5,1) =   hm1**2 / 2.d0
  m(5,2) =   0.d0
  m(5,3) =   hp1**2 / 2.d0
  m(5,4) = - hm2**4 / 24.d0
  m(5,5) = - hm1**4 / 24.d0
  m(5,6) =   0.d0
  m(5,7) = - hp1**4 / 24.d0
  m(5,8) = - hp2**4 / 24.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - 1.d0
  m(6,5) = - 1.d0
  m(6,6) = - 1.d0
  m(6,7) = - 1.d0
  m(6,8) = - 1.d0
 
  m(7,1) = - hm1**3 / 6.d0
  m(7,2) =   0.d0
  m(7,3) =   hp1**3 / 6.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) = - hp1**5 / 120.d0
  m(7,8) = - hp2**5 / 120.d0

  m(8,1) =   hm1**4 / 24.d0
  m(8,2) =   0.d0
  m(8,3) =   hp1**4 / 24.d0
  m(8,4) = - hm2**6 / 720.d0
  m(8,5) = - hm1**6 / 720.d0
  m(8,6) =   0.d0
  m(8,7) = - hp1**6 / 720.d0
  m(8,8) = - hp2**6 / 720.d0
 
  return

 end subroutine center_point_second_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_center_point_second_derivative_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(2) = 1.d0
 
  return

 end subroutine rhs_center_point_second_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine penultimate_point_second_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=8), intent(out)  :: m(n,n)
  real(kind=8)               :: hp1, hm1, hm2, hm3, hm4, hm5, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 =       stf**5
  hm1 =       stf**4
  hm2 = hm1 + stf**3
  hm3 = hm2 + stf**2
  hm4 = hm3 + stf
  hm5 = hm4 + dy

  ! matrix definition
  m(1,1)  =   1.d0
  m(1,2)  =   0.d0
  m(1,3)  =   0.d0
  m(1,4)  =   0.d0
  m(1,5)  =   0.d0
  m(1,6)  =   0.d0
  m(1,7)  =   0.d0
  m(1,8)  =   0.d0
  m(1,9)  =   0.d0
  m(1,10) =   0.d0
 
  m(2,1)  =   1.d0
  m(2,2)  =   1.d0
  m(2,3)  =   1.d0
  m(2,4)  = - hp1**2 / 2.d0
  m(2,5)  =   0.d0
  m(2,6)  = - hm1**2 / 2.d0
  m(2,7)  = - hm2**2 / 2.d0
  m(2,8)  = - hm3**2 / 2.d0
  m(2,9)  = - hm4**2 / 2.d0
  m(2,10) = - hm5**2 / 2.d0
 
  m(3,1)  = - hp1
  m(3,2)  =   0.d0
  m(3,3)  =   hm1
  m(3,4)  =   hp1**3 / 6.d0
  m(3,5)  =   0.d0
  m(3,6)  = - hm1**3 / 6.d0
  m(3,7)  = - hm2**3 / 6.d0
  m(3,8)  = - hm3**3 / 6.d0
  m(3,9)  = - hm4**3 / 6.d0
  m(3,10) = - hm5**3 / 6.d0
 
  m(4,1)  =   hp1**2 / 2.d0
  m(4,2)  =   0.d0
  m(4,3)  =   hm1**2 / 2.d0
  m(4,4)  = - hp1**4 / 24.d0
  m(4,5)  =   0.d0
  m(4,6)  = - hm1**4 / 24.d0
  m(4,7)  = - hm2**4 / 24.d0
  m(4,8)  = - hm3**4 / 24.d0
  m(4,9)  = - hm4**4 / 24.d0
  m(4,10) = - hm5**4 / 24.d0
 
  m(5,1)  =   0.d0
  m(5,2)  =   0.d0
  m(5,3)  =   0.d0
  m(5,4)  = - 1.d0
  m(5,5)  = - 1.d0
  m(5,6)  = - 1.d0
  m(5,7)  = - 1.d0
  m(5,8)  = - 1.d0
  m(5,9)  = - 1.d0
  m(5,10) = - 1.d0

  m(6,1)  =   0.d0
  m(6,2)  =   0.d0
  m(6,3)  =   0.d0
  m(6,4)  =   0.d0
  m(6,5)  =   0.d0
  m(6,6)  = - 1.d0
  m(6,7)  =   0.d0
  m(6,8)  =   0.d0
  m(6,9)  =   0.d0
  m(6,10) =   0.d0
 
  m(7,1)  = - hp1**3 / 6.d0
  m(7,2)  =   0.d0
  m(7,3)  =   hm1**3 / 6.d0
  m(7,4)  =   hp1**5 / 120.d0
  m(7,5)  =   0.d0
  m(7,6)  = - hm1**5 / 120.d0
  m(7,7)  = - hm2**5 / 120.d0
  m(7,8)  = - hm3**5 / 120.d0
  m(7,9)  = - hm4**5 / 120.d0
  m(7,10) = - hm5**5 / 120.d0
 
  m(8,1)  =   hp1**4 / 24.d0
  m(8,2)  =   0.d0
  m(8,3)  =   hm1**4 / 24.d0
  m(8,4)  = - hp1**6 / 720.d0
  m(8,5)  =   0.d0
  m(8,6)  = - hm1**6 / 720.d0
  m(8,7)  = - hm2**6 / 720.d0
  m(8,8)  = - hm3**6 / 720.d0
  m(8,9)  = - hm4**6 / 720.d0
  m(8,10) = - hm5**6 / 720.d0
 
  m(9,1)  = - hp1**5 / 120.d0
  m(9,2)  =   0.d0
  m(9,3)  =   hm1**5 / 120.d0
  m(9,4)  =   hp1**7 / 5040.d0
  m(9,5)  =   0.d0
  m(9,6)  = - hm1**7 / 5040.d0
  m(9,7)  = - hm2**7 / 5040.d0
  m(9,8)  = - hm3**7 / 5040.d0
  m(9,9)  = - hm4**7 / 5040.d0
  m(9,10) = - hm5**7 / 5040.d0
 
  m(10,1)  =   0.d0
  m(10,2)  =   0.d0
  m(10,3)  =   0.d0
  m(10,4)  =   hp1
  m(10,5)  =   0.d0
  m(10,6)  = - hm1
  m(10,7)  = - hm2
  m(10,8)  = - hm3
  m(10,9)  = - hm4
  m(10,10) = - hm5

  return

 end subroutine penultimate_point_second_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_penultimate_point_second_derivative_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 2, n
   b(i) = 0.d0
  end do
  b(1) =   1.d0
  b(6) = - 5.25d0
 
  return

  end subroutine rhs_penultimate_point_second_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine last_point_second_derivative_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm1, hm2, hm3, hm4, hm5, dy

  dy = 1.d0
  ! 'normalized' length
  hm1 =       stf**4
  hm2 = hm1 + stf**3
  hm3 = hm2 + stf**2
  hm4 = hm3 + stf
  hm5 = hm4 + dy

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  m(1,8) =   0.d0
 
  m(2,1) =   1.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) = - hm1**2 / 2.d0
  m(2,5) = - hm2**2 / 2.d0
  m(2,6) = - hm3**2 / 2.d0
  m(2,7) = - hm4**2 / 2.d0
  m(2,8) = - hm5**2 / 2.d0
 
  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
  m(3,5) = - 1.d0
  m(3,6) = - 1.d0
  m(3,7) = - 1.d0
  m(3,8) = - 1.d0

  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - hm1
  m(6,5) = - hm2
  m(6,6) = - hm3
  m(6,7) = - hm4
  m(6,8) = - hm5
 
  m(4,1) =   0.d0
  m(4,2) =   hm1
  m(4,3) =   0.d0
  m(4,4) = - hm1**3 / 6.d0
  m(4,5) = - hm2**3 / 6.d0
  m(4,6) = - hm3**3 / 6.d0
  m(4,7) = - hm4**3 / 6.d0
  m(4,8) = - hm5**3 / 6.d0

  m(5,1) =   0.d0
  m(5,2) =   hm1**2 / 2.d0
  m(5,3) =   0.d0
  m(5,4) = - hm1**4 / 24.d0
  m(5,5) = - hm2**4 / 24.d0
  m(5,6) = - hm3**4 / 24.d0
  m(5,7) = - hm4**4 / 24.d0
  m(5,8) = - hm5**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hm1**3 / 6.d0
  m(7,3) =   0.d0
  m(7,4) = - hm1**5 / 120.d0
  m(7,5) = - hm2**5 / 120.d0
  m(7,6) = - hm3**5 / 120.d0
  m(7,7) = - hm4**5 / 120.d0
  m(7,8) = - hm5**5 / 120.d0

  m(8,1) =   0.d0
  m(8,2) =   hm1**4 / 24.d0
  m(8,3) =   0.d0
  m(8,4) = - hm1**6 / 720.d0
  m(8,5) = - hm2**6 / 720.d0
  m(8,6) = - hm3**6 / 720.d0
  m(8,7) = - hm4**6 / 720.d0
  m(8,8) = - hm5**6 / 720.d0

  return

 end subroutine last_point_second_derivative_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_last_point_second_derivative_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)

  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0

  return

  end subroutine rhs_last_point_second_derivative_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine second_point_poisson_matrix_coef(m, strl, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=8), intent(in)   :: strl
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: hm1, hp1, hp2, dy

  dy = 1.d0
  ! 'normalized' length
  hm1 = dy
  hp1 =       strl    * dy
  hp2 = hp1 + strl**2 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   hm1
  m(1,5) =   0.d0
  m(1,6) = - hp1
  m(1,7) = - hp2
 
  m(2,1) =   0.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) =   0.d0
  m(2,5) =   0.d0
  m(2,6) =   0.d0
  m(2,7) =   0.d0
 
  m(3,1) = - hm1
  m(3,2) =   1.d0
  m(3,3) =   1.d0
  m(3,4) = - hm1**2 / 2.d0
  m(3,5) =   0.d0
  m(3,6) = - hp1**2 / 2.d0
  m(3,7) = - hp2**2 / 2.d0
 
  m(4,1) =   hm1**2 / 2.d0
  m(4,2) =   0.d0
  m(4,3) =   hp1
  m(4,4) =   hm1**3 / 6.d0
  m(4,5) =   0.d0
  m(4,6) = - hp1**3 / 6.d0
  m(4,7) = - hp2**3 / 6.d0
 
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) =   0.d0
  m(5,4) = - 1.d0
  m(5,5) = - 1.d0
  m(5,6) = - 1.d0
  m(5,7) = - 1.d0
 
  m(6,1) = - hm1**3 / 6.d0
  m(6,2) =   0.d0
  m(6,3) =   hp1**2 / 2.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hp1**4 / 24.d0
  m(6,7) = - hp2**4 / 24.d0

  m(7,1) =   hm1**4 / 24.d0
  m(7,2) =   0.d0
  m(7,3) =   hp1**3 / 6.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) = - hp1**5 / 120.d0
  m(7,7) = - hp2**5 / 120.d0

 return

 end subroutine second_point_poisson_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_second_point_poisson_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=16), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(2) = 4.d0 / 5.d0
 
  return

 end subroutine rhs_second_point_poisson_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine center_point_poisson_matrix_coef(m, strl, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=8), intent(in)   :: strl
  real(kind=16)              :: hp1, hp2, hm1, hm2
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: dy

  dy = 1.d0
  ! 'normalized' length
  hm1 =       strl
  hm2 = dy  + strl
  hp1 =       strl**2
  hp2 = hp1 + strl**3

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   1.d0
  m(1,3) =   1.d0
  m(1,4) = - hm2**2 / 2.d0
  m(1,5) = - hm1**2 / 2.d0
  m(1,6) =   0.d0
  m(1,7) = - hp1**2 / 2.d0
  m(1,8) = - hp2**2 / 2.d0
 
  m(2,1) =   0.d0
  m(2,2) =   1.d0
  m(2,3) =   0.d0
  m(2,4) =   0.d0
  m(2,5) =   0.d0
  m(2,6) =   0.d0
  m(2,7) =   0.d0
  m(2,8) =   0.d0
 
  m(3,1) = - hm1
  m(3,2) =   0.d0
  m(3,3) =   hp1
  m(3,4) =   hm2**3 / 6.d0
  m(3,5) =   hm1**3 / 6.d0
  m(3,6) =   0.d0
  m(3,7) = - hp1**3 / 6.d0
  m(3,8) = - hp2**3 / 6.d0
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) =   hm2
  m(4,5) =   hm1
  m(4,6) =   0.d0
  m(4,7) = - hp1
  m(4,8) = - hp2
 
  m(5,1) =   hm1**2 / 2.d0
  m(5,2) =   0.d0
  m(5,3) =   hp1**2 / 2.d0
  m(5,4) = - hm2**4 / 24.d0
  m(5,5) = - hm1**4 / 24.d0
  m(5,6) =   0.d0
  m(5,7) = - hp1**4 / 24.d0
  m(5,8) = - hp2**4 / 24.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - 1.d0
  m(6,5) = - 1.d0
  m(6,6) = - 1.d0
  m(6,7) = - 1.d0
  m(6,8) = - 1.d0
 
  m(7,1) = - hm1**3 / 6.d0
  m(7,2) =   0.d0
  m(7,3) =   hp1**3 / 6.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) = - hp1**5 / 120.d0
  m(7,8) = - hp2**5 / 120.d0

  m(8,1) =   hm1**4 / 24.d0
  m(8,2) =   0.d0
  m(8,3) =   hp1**4 / 24.d0
  m(8,4) = - hm2**6 / 720.d0
  m(8,5) = - hm1**6 / 720.d0
  m(8,6) =   0.d0
  m(8,7) = - hp1**6 / 720.d0
  m(8,8) = - hp2**6 / 720.d0
 
  return

 end subroutine center_point_poisson_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_center_point_poisson_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=16), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(2) = 11.d0/ 15.d0
 
  return

 end subroutine rhs_center_point_poisson_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine penultimate_point_poisson_matrix_coef(m, strl, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=8), intent(in)   :: strl
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: hp1, hm1, dy

  dy = 1.d0
  ! 'normalized' length
  hp1 = strl
  hm1 = dy
 
  ! matrix definition
  m(4,1) =   1.d0
  m(4,2) = - hm1**2 / 2.d0
  m(4,3) =   0.d0
  m(4,4) = - hp1**2 / 2.d0
 
  m(2,1) =   0.d0
  m(2,2) =   hm1
  m(2,3) =   0.d0
  m(2,4) = - hp1
 
  m(3,1) =   0.d0
  m(3,2) = - 1.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
 
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
 
  return

 end subroutine penultimate_point_poisson_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_penultimate_point_poisson_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)        :: n
  integer                    :: i
  real(kind=16), intent(out) :: b(n)
 
  do i = 2, n
   b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

  end subroutine rhs_penultimate_point_poisson_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine last_point_poisson_matrix_coef(m, strl, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=8), intent(in)   :: strl
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: hm1, hm2, dy

  dy = 1.d0
  ! 'normalized' length
  hm1 = strl
  hm2 = strl + dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) =   0.d0
  m(2,4) =   hm1
  m(2,5) =   hm2
 
  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) = - 1.d0
  m(3,4) = - 1.d0
  m(3,5) = - 1.d0

  m(4,1) =   1.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) = - hm1**2 / 2.d0
  m(4,5) = - hm2**2 / 2.d0

  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) =   0.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   hm2**3 / 6.d0

  return

 end subroutine last_point_poisson_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_last_point_poisson_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)        :: n
  integer                    :: i
  real(kind=16), intent(out) :: b(n)

  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0

  return

  end subroutine rhs_last_point_poisson_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine vorticity_at_wall_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: hp1, hp2, hp3, hp4, hp5, hp6, dy

  dy = 1.d0
  ! 'normalized' length
  hp1 = dy
  hp2 = hp1 + stf    * dy
  hp3 = hp2 + stf**2 * dy
  hp4 = hp3 + stf**3 * dy
  hp5 = hp4 + stf**4 * dy
  hp6 = hp5 + stf**5 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) = - hp1**2 / 2.d0
  m(1,5) = - hp2**2 / 2.d0
  m(1,6) = - hp3**2 / 2.d0
  m(1,7) = - hp4**2 / 2.d0
  m(1,8) = - hp5**2 / 2.d0
  m(1,9) = - hp6**2 / 2.d0

  m(2,1) =   0.d0
  m(2,2) =   1.d0
  m(2,3) = - 1.d0
  m(2,4) = - hp1
  m(2,5) = - hp2
  m(2,6) = - hp3
  m(2,7) = - hp4
  m(2,8) = - hp5
  m(2,9) = - hp6
 
  m(3,1) =   0.d0
  m(3,2) =   0.d0
  m(3,3) =   1.d0
  m(3,4) =   1.d0
  m(3,5) =   1.d0
  m(3,6) =   1.d0
  m(3,7) =   1.d0
  m(3,8) =   1.d0
  m(3,9) =   1.d0

  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) =   0.d0
  m(4,4) = - hp1**3 / 6.d0
  m(4,5) = - hp2**3 / 6.d0
  m(4,6) = - hp3**3 / 6.d0
  m(4,7) = - hp4**3 / 6.d0
  m(4,8) = - hp5**3 / 6.d0
  m(4,9) = - hp6**3 / 6.d0
 
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) =   0.d0
  m(5,4) = - hp1**4 / 24.d0
  m(5,5) = - hp2**4 / 24.d0
  m(5,6) = - hp3**4 / 24.d0
  m(5,7) = - hp4**4 / 24.d0
  m(5,8) = - hp5**4 / 24.d0
  m(5,9) = - hp6**4 / 24.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) =   0.d0
  m(6,4) = - hp1**5 / 120.d0
  m(6,5) = - hp2**5 / 120.d0
  m(6,6) = - hp3**5 / 120.d0
  m(6,7) = - hp4**5 / 120.d0
  m(6,8) = - hp5**5 / 120.d0
  m(6,9) = - hp6**5 / 120.d0
 
  m(7,1) =   0.d0
  m(7,2) =   0.d0
  m(7,3) =   0.d0
  m(7,4) = - hp1**6 / 720.d0
  m(7,5) = - hp2**6 / 720.d0
  m(7,6) = - hp3**6 / 720.d0
  m(7,7) = - hp4**6 / 720.d0
  m(7,8) = - hp5**6 / 720.d0
  m(7,9) = - hp6**6 / 720.d0
 
  m(8,1) =   0.d0
  m(8,2) =   0.d0
  m(8,3) =   0.d0
  m(8,4) = - hp1**7 / 5040.d0
  m(8,5) = - hp2**7 / 5040.d0
  m(8,6) = - hp3**7 / 5040.d0
  m(8,7) = - hp4**7 / 5040.d0
  m(8,8) = - hp5**7 / 5040.d0
  m(8,9) = - hp6**7 / 5040.d0
 
! m(9,1) =   0.d0
! m(9,2) =   0.d0
! m(9,3) =   0.d0
! m(9,4) = - hp1**8 / 40320.d0
! m(9,5) = - hp2**8 / 40320.d0
! m(9,6) = - hp3**8 / 40320.d0
! m(9,7) = - hp4**8 / 40320.d0
! m(9,8) = - hp5**8 / 40320.d0
! m(9,9) = - hp6**8 / 40320.d0

  m(9,1) =   0.d0
  m(9,2) =   0.d0
  m(9,3) =   0.d0
  m(9,4) =   0.d0
  m(9,5) =   0.d0
  m(9,6) =   0.d0
  m(9,7) =   0.d0
  m(9,8) =   0.d0
  m(9,9) =   1.d0

 return

 end subroutine vorticity_at_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_vorticity_at_wall_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)        :: n
  integer                    :: i
  real(kind=16), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(9) = - 1.d0 / 18.d0
 
  return
 
 end subroutine rhs_vorticity_at_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine derivative_coefficient_generation
 
  implicit none
  real(kind=8) :: fp_fd_c(7,7), fp_fd_b(7)
  real(kind=8) :: fp_fd_c_e(6,6), fp_fd_b_e(6)
  real(kind=8) :: sp_fd_c(9,9), sp_fd_b(9)
  real(kind=8) :: cp_fd_c(8,8), cp_fd_b(8)
  real(kind=8) :: pp_fd_c(9,9), pp_fd_b(9)
  real(kind=8) :: lp_fd_c(7,7), lp_fd_b(7)

  real(kind=8) :: fp_sd_c(8,8)  , fp_sd_b(8)
  real(kind=8) :: sp_sd_c(10,10), sp_sd_b(10)
  real(kind=8) :: cp_sd_c(8,8)  , cp_sd_b(8)
  real(kind=8) :: pp_sd_c(10,10), pp_sd_b(10)
  real(kind=8) :: lp_sd_c(8,8)  , lp_sd_b(8)
 
  real(kind=8) :: sp_poi_c(7,7), sp_poi_b(7)
  real(kind=8) :: cp_poi_c(8,8), cp_poi_b(8)
  real(kind=8) :: pp_poi_c(4,4), pp_poi_b(4)
  real(kind=8) :: lp_poi_c(5,5), lp_poi_b(5)
 
  real(kind=16) :: w_at_w_c(9,9), w_at_w_b(9), w_at_w_r(9)
 
  real(kind=16) :: dwydy_c(8,8), dwydy_b(8), dwydy_r(8)

  ! First derivative coefficients
  call first_point_first_derivative_matrix_coef(fp_fd_c, 7)
  call rhs_first_point_first_derivative_matrix_coef(fp_fd_b, 7)
  call ludecomp(fp_fd_c, fp_fd_b, fp_fd_coef, 7)

  call first_point_first_derivative_matrix_coef_explicit(fp_fd_c_e, 6)
  call rhs_first_point_first_derivative_matrix_coef_explicit(fp_fd_b_e, 6)
  call ludecomp(fp_fd_c_e, fp_fd_b_e, fp_fd_coef_e, 6)

  call second_point_first_derivative_matrix_coef(sp_fd_c, 9)
  call rhs_second_point_first_derivative_matrix_coef(sp_fd_b, 9)
  call ludecomp(sp_fd_c, sp_fd_b, sp_fd_coef, 9)
 
  call center_point_first_derivative_matrix_coef(cp_fd_c, 8)
  call rhs_center_point_first_derivative_matrix_coef(cp_fd_b, 8)
  call ludecomp(cp_fd_c, cp_fd_b, cp_fd_coef, 8)
 
  call penultimate_point_first_derivative_matrix_coef(pp_fd_c, 9)
  call rhs_penultimate_point_first_derivative_matrix_coef(pp_fd_b, 9)
  call ludecomp(pp_fd_c, pp_fd_b, pp_fd_coef, 9)
 
  call last_point_first_derivative_matrix_coef(lp_fd_c, 7)
  call rhs_last_point_first_derivative_matrix_coef(lp_fd_b, 7)
  call ludecomp(lp_fd_c, lp_fd_b, lp_fd_coef, 7)
 
  ! Second derivative coefficients
  call first_point_second_derivative_matrix_coef(fp_sd_c, 8)
  call rhs_first_point_second_derivative_matrix_coef(fp_sd_b, 8)
  call ludecomp(fp_sd_c, fp_sd_b, fp_sd_coef, 8)
 
  call second_point_second_derivative_matrix_coef(sp_sd_c, 10)
  call rhs_second_point_second_derivative_matrix_coef(sp_sd_b, 10)
  call ludecomp(sp_sd_c, sp_sd_b, sp_sd_coef, 10)
 
  call center_point_second_derivative_matrix_coef(cp_sd_c, 8)
  call rhs_center_point_second_derivative_matrix_coef(cp_sd_b, 8)
  call ludecomp(cp_sd_c, cp_sd_b, cp_sd_coef, 8)
 
  call penultimate_point_second_derivative_matrix_coef(pp_sd_c, 10)
  call rhs_penultimate_point_second_derivative_matrix_coef(pp_sd_b, 10)
  call ludecomp(pp_sd_c, pp_sd_b, pp_sd_coef, 10)
 
  call last_point_second_derivative_matrix_coef(lp_sd_c, 8)
  call rhs_last_point_second_derivative_matrix_coef(lp_sd_b, 8)
  call ludecomp(lp_sd_c, lp_sd_b, lp_sd_coef, 8)
 
  call vorticity_at_wall_matrix_coef(w_at_w_c, 9)
  call rhs_vorticity_at_wall_matrix_coef(w_at_w_b, 9)
  call ludecomp_hp(w_at_w_c, w_at_w_b, w_at_w_r, 9)
 
  call dwydy_matrix_coef(dwydy_c, 8)
  call rhs_dwydy_matrix_coef(dwydy_b, 8)
  call ludecomp_hp(dwydy_c, dwydy_b, dwydy_r, 8)
  w_at_w_coef = w_at_w_r ! vector and casting operation
  dwydy_coef  = dwydy_r ! vector and casting operation

 end subroutine derivative_coefficient_generation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine poisson_coefficient_generation(stf_lvl)

  implicit none
  real(kind=8), intent(in) :: stf_lvl

  real(kind=16) :: sp_poi_c(7,7), sp_poi_b(7), sp_poi_r(7)
  real(kind=16) :: cp_poi_c(8,8), cp_poi_b(8), cp_poi_r(8)
  real(kind=16) :: pp_poi_c(4,4), pp_poi_b(4), pp_poi_r(4)
  real(kind=16) :: lp_poi_c(5,5), lp_poi_b(5), lp_poi_r(5)
 
  ! Poisson derivative coefficients
  call second_point_poisson_matrix_coef(sp_poi_c, stf_lvl, 7)
  call rhs_second_point_poisson_matrix_coef(sp_poi_b, 7)
  call ludecomp_hp(sp_poi_c, sp_poi_b, sp_poi_r, 7)
 
  call center_point_poisson_matrix_coef(cp_poi_c, stf_lvl, 8)
  call rhs_center_point_poisson_matrix_coef(cp_poi_b, 8)
  call ludecomp_hp(cp_poi_c, cp_poi_b, cp_poi_r, 8)
 
  call penultimate_point_poisson_matrix_coef(pp_poi_c, stf_lvl, 4)
  call rhs_penultimate_point_poisson_matrix_coef(pp_poi_b, 4)
  call ludecomp_hp(pp_poi_c, pp_poi_b, pp_poi_r, 4)
 
  call last_point_poisson_matrix_coef(lp_poi_c, stf_lvl, 5)
  call rhs_last_point_poisson_matrix_coef(lp_poi_b, 5)
  call ludecomp_hp(lp_poi_c, lp_poi_b, lp_poi_r, 5)
 
  sp_poi_coef = sp_poi_r
  cp_poi_coef = cp_poi_r
  lp_poi_coef = lp_poi_r
  pp_poi_coef = pp_poi_r
 
 end subroutine poisson_coefficient_generation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine dwydy_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)        :: n
  real(kind=16), intent(out) :: m(n,n)
  real(kind=16)              :: hp1, hp2, hp3, hp4, hp5, hp6, dy

  dy = 1.d0
  ! 'normalized' length
  hp1 = dy
  hp2 = hp1 + stf**1
  hp3 = hp2 + stf**2
  hp4 = hp3 + stf**3
  hp5 = hp4 + stf**4
  hp6 = hp5 + stf**5
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  m(1,8) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  m(2,8) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) =   0.d0
  m(3,3) = - hp1
  m(3,4) = - hp2
  m(3,5) = - hp3
  m(3,6) = - hp4
  m(3,7) = - hp5
  m(3,8) = - hp6
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) = - hp1**2 / 2.d0
  m(4,4) = - hp2**2 / 2.d0
  m(4,5) = - hp3**2 / 2.d0
  m(4,6) = - hp4**2 / 2.d0
  m(4,7) = - hp5**2 / 2.d0
  m(4,8) = - hp6**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) = - hp1**3 / 6.d0
  m(5,4) = - hp2**3 / 6.d0
  m(5,5) = - hp3**3 / 6.d0
  m(5,6) = - hp4**3 / 6.d0
  m(5,7) = - hp5**3 / 6.d0
  m(5,8) = - hp6**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) = - hp1**4 / 24.d0
  m(6,4) = - hp2**4 / 24.d0
  m(6,5) = - hp3**4 / 24.d0
  m(6,6) = - hp4**4 / 24.d0
  m(6,7) = - hp5**4 / 24.d0
  m(6,8) = - hp6**4 / 24.d0

  m(7,1) =   0.d0
  m(7,2) =   0.d0
  m(7,3) = - hp1**5 / 120.d0
  m(7,4) = - hp2**5 / 120.d0
  m(7,5) = - hp3**5 / 120.d0
  m(7,6) = - hp4**5 / 120.d0
  m(7,7) = - hp5**5 / 120.d0
  m(7,8) = - hp6**5 / 120.d0
 
  m(8,1) =   0.d0
  m(8,2) =   0.d0
  m(8,3) = - hp1**6 / 720.d0
  m(8,4) = - hp2**6 / 720.d0
  m(8,5) = - hp3**6 / 720.d0
  m(8,6) = - hp4**6 / 720.d0
  m(8,7) = - hp5**6 / 720.d0
  m(8,8) = - hp6**6 / 720.d0

 return

 end subroutine dwydy_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine rhs_dwydy_matrix_coef(b, n)

  ! right-hand side of second derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=16), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0

  return

 end subroutine rhs_dwydy_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 subroutine first_point_first_derivative_matrix_coef_explicit(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hp3, hp4, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 = dy
  hp2 = hp1 + stf    * dy
  hp3 = hp2 + stf**2 * dy
  hp4 = hp3 + stf**3 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0

  m(3,1) =   1.d0
  m(3,2) =   0.d0
  m(3,3) = - hp1
  m(3,4) = - hp2
  m(3,5) = - hp3
  m(3,6) = - hp4

  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) = - hp1**2 / 2.d0
  m(4,4) = - hp2**2 / 2.d0
  m(4,5) = - hp3**2 / 2.d0
  m(4,6) = - hp4**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) = - hp1**3 / 6.d0
  m(5,4) = - hp2**3 / 6.d0
  m(5,5) = - hp3**3 / 6.d0
  m(5,6) = - hp4**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) = - hp1**4 / 24.d0
  m(6,4) = - hp2**4 / 24.d0
  m(6,5) = - hp3**4 / 24.d0
  m(6,6) = - hp4**4 / 24.d0

 return

 end subroutine first_point_first_derivative_matrix_coef_explicit
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_first_point_first_derivative_matrix_coef_explicit(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 2, n
    b(i) = 0.d0
  end do
  b(1) = 12.d0
 
  return

 end subroutine rhs_first_point_first_derivative_matrix_coef_explicit

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

end module derivative

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
