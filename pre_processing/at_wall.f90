
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

module atwall_calculation
        
 use derivative_commom_variables
 use solvers

 contains

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine first_point_wall_matrix_coef(m, n)
  
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
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   0.d0
  m(3,3) = - hp1
  m(3,4) = - hp2
  m(3,5) = - hp3
  m(3,6) = - hp4
  m(3,7) = - hp5
  
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) = - hp1**2 / 2.d0
  m(4,4) = - hp2**2 / 2.d0
  m(4,5) = - hp3**2 / 2.d0
  m(4,6) = - hp4**2 / 2.d0
  m(4,7) = - hp5**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) = - hp1**3 / 6.d0
  m(5,4) = - hp2**3 / 6.d0
  m(5,5) = - hp3**3 / 6.d0
  m(5,6) = - hp4**3 / 6.d0
  m(5,7) = - hp5**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) = - hp1**4 / 24.d0
  m(6,4) = - hp2**4 / 24.d0
  m(6,5) = - hp3**4 / 24.d0
  m(6,6) = - hp4**4 / 24.d0
  m(6,7) = - hp5**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   0.d0
  m(7,3) = - hp1**5 / 120.d0
  m(7,4) = - hp2**5 / 120.d0
  m(7,5) = - hp3**5 / 120.d0
  m(7,6) = - hp4**5 / 120.d0
  m(7,7) = - hp5**5 / 120.d0

 return

 end subroutine first_point_wall_matrix_coef
   
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
     
 subroutine rhs_first_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_first_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine second_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm1, hp1, hp2, hp3, hp4, dy
  
  dy = 1.d0
  ! 'normalized' length
  hm1 = dy
  hp1 = stf * dy
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
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   hm1
  m(3,3) =   0.d0
  m(3,4) = - hp1
  m(3,5) = - hp2
  m(3,6) = - hp3
  m(3,7) = - hp4
  
  m(4,1) = 0.d0
  m(4,2) = -hm1**2 / 2.d0
  m(4,3) =  0.d0
  m(4,4) = -hp1**2 / 2.d0
  m(4,5) = -hp2**2 / 2.d0
  m(4,6) = -hp3**2 / 2.d0
  m(4,7) = -hp4**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   hm1**3 / 6.d0
  m(5,3) =   0.d0
  m(5,4) = - hp1**3 / 6.d0
  m(5,5) = - hp2**3 / 6.d0
  m(5,6) = - hp3**3 / 6.d0
  m(5,7) = - hp4**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) = - hm1**4 / 24.d0
  m(6,3) =   0.d0
  m(6,4) = - hp1**4 / 24.d0
  m(6,5) = - hp2**4 / 24.d0
  m(6,6) = - hp3**4 / 24.d0
  m(6,7) = - hp4**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   hm1**5 / 120.d0
  m(7,3) =   0.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) = - hp2**5 / 120.d0
  m(7,6) = - hp3**5 / 120.d0
  m(7,7) = - hp4**5 / 120.d0

 return

 end subroutine second_point_wall_matrix_coef
   
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
     
 subroutine rhs_second_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_second_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 subroutine third_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm2, hm1, hp1, hp2, hp3, dy
  
  dy = 1.d0
  ! 'normalized' length
  hm2 = dy  + stf    * dy
  hm1 =       stf    * dy
  hp1 =       stf**2 * dy
  hp2 = hp1 + stf**3 * dy
  hp3 = hp2 + stf**4 * dy
  
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   hm2
  m(3,3) =   hm1
  m(3,4) =   0.d0
  m(3,5) = - hp1
  m(3,6) = - hp2
  m(3,7) = - hp3
  
  m(4,1) =   0.d0
  m(4,2) = - hm2**2 / 2.d0
  m(4,3) = - hm1**2 / 2.d0
  m(4,4) =   0.d0
  m(4,5) = - hp1**2 / 2.d0
  m(4,6) = - hp2**2 / 2.d0
  m(4,7) = - hp3**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   hm2**3 / 6.d0
  m(5,3) =   hm1**3 / 6.d0
  m(5,4) =   0.d0
  m(5,5) = - hp1**3 / 6.d0
  m(5,6) = - hp2**3 / 6.d0
  m(5,7) = - hp3**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) = - hm2**4 / 24.d0
  m(6,3) = - hm1**4 / 24.d0
  m(6,4) =   0.d0
  m(6,5) = - hp1**4 / 24.d0
  m(6,6) = - hp2**4 / 24.d0
  m(6,7) = - hp3**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   hm2**5 / 120.d0
  m(7,3) =   hm1**5 / 120.d0
  m(7,4) =   0.d0
  m(7,5) = - hp1**5 / 120.d0
  m(7,6) = - hp2**5 / 120.d0
  m(7,7) = - hp3**5 / 120.d0

 return

 end subroutine third_point_wall_matrix_coef
   
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_third_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_third_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine fourth_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm3, hm2, hm1, hp1, hp2, dy
  
  dy = 1.d0
  ! 'normalized' length
  hm1 =       stf**2 * dy
  hm2 = hm1 + stf    * dy
  hm3 = hm2 + dy
  hp1 =       stf**3 * dy
  hp2 = hp1 + stf**4 * dy
  
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   hm3
  m(3,3) =   hm2
  m(3,4) =   hm1
  m(3,5) =   0.d0
  m(3,6) = - hp1
  m(3,7) = - hp2
  
  m(4,1) =   0.d0
  m(4,2) = - hm3**2 / 2.d0
  m(4,3) = - hm2**2 / 2.d0
  m(4,4) = - hm1**2 / 2.d0
  m(4,5) =   0.d0
  m(4,6) = - hp1**2 / 2.d0
  m(4,7) = - hp2**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   hm3**3 / 6.d0
  m(5,3) =   hm2**3 / 6.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   0.d0
  m(5,6) = - hp1**3 / 6.d0
  m(5,7) = - hp2**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) = - hm3**4 / 24.d0
  m(6,3) = - hm2**4 / 24.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hp1**4 / 24.d0
  m(6,7) = - hp2**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   hm3**5 / 120.d0
  m(7,3) =   hm2**5 / 120.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) = - hp1**5 / 120.d0
  m(7,7) = - hp2**5 / 120.d0

 return

 end subroutine fourth_point_wall_matrix_coef
   
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_fourth_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_fourth_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine fifth_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm4, hm3, hm2, hm1, hp1, dy
  
  dy = 1.d0
  ! 'normalized' length
  hm1 =       stf**3 * dy
  hm2 = hm1 + stf**2 * dy
  hm3 = hm2 + stf    * dy
  hm4 = hm3 + dy
  hp1 =       stf**4 * dy
  
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   hm4
  m(3,3) =   hm3
  m(3,4) =   hm2
  m(3,5) =   hm1
  m(3,6) =   0.d0
  m(3,7) = - hp1
  
  m(4,1) =   0.d0
  m(4,2) = - hm4**2 / 2.d0
  m(4,3) = - hm3**2 / 2.d0
  m(4,4) = - hm2**2 / 2.d0
  m(4,5) = - hm1**2 / 2.d0
  m(4,6) =   0.d0
  m(4,7) = - hp1**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   hm4**3 / 6.d0
  m(5,3) =   hm3**3 / 6.d0
  m(5,4) =   hm2**3 / 6.d0
  m(5,5) =   hm1**3 / 6.d0
  m(5,6) =   0.d0
  m(5,7) = - hp1**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) = - hm4**4 / 24.d0
  m(6,3) = - hm3**4 / 24.d0
  m(6,4) = - hm2**4 / 24.d0
  m(6,5) = - hm1**4 / 24.d0
  m(6,6) =   0.d0
  m(6,7) = - hp1**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   hm4**5 / 120.d0
  m(7,3) =   hm3**5 / 120.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) = - hp1**5 / 120.d0

 return

 end subroutine fifth_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_fifth_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_fifth_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine jmax_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm5, hm4, hm3, hm2, hm1, dy
  
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
  
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
  
  m(3,1) =   1.d0
  m(3,2) =   0.d0
  m(3,3) =   hm1
  m(3,4) =   hm2
  m(3,5) =   hm3
  m(3,6) =   hm4
  m(3,7) =   hm5
  
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) = - hm1**2 / 2.d0
  m(4,4) = - hm2**2 / 2.d0
  m(4,5) = - hm3**2 / 2.d0
  m(4,6) = - hm4**2 / 2.d0
  m(4,7) = - hm5**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) =   hm1**3 / 6.d0
  m(5,4) =   hm2**3 / 6.d0
  m(5,5) =   hm3**3 / 6.d0
  m(5,6) =   hm4**3 / 6.d0
  m(5,7) =   hm5**3 / 6.d0
  
  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) = - hm1**4 / 24.d0
  m(6,4) = - hm2**4 / 24.d0
  m(6,5) = - hm3**4 / 24.d0
  m(6,6) = - hm4**4 / 24.d0
  m(6,7) = - hm5**4 / 24.d0
  
  m(7,1) =   0.d0
  m(7,2) =   0.d0
  m(7,3) =   hm1**5 / 120.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm3**5 / 120.d0
  m(7,6) =   hm4**5 / 120.d0
  m(7,7) =   hm5**5 / 120.d0

 return

 end subroutine jmax_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_jmax_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_jmax_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine jmaxm1_point_wall_matrix_coef(m, n)

  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hm4, hm3, hm2, hm1, dy
  
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

  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) = - hp1
  m(3,3) =   0.d0
  m(3,4) =   hm1
  m(3,5) =   hm2
  m(3,6) =   hm3
  m(3,7) =   hm4
 
  m(4,1) =   0.d0
  m(4,2) = - hp1**2 / 2.d0
  m(4,3) =   0.d0
  m(4,4) = - hm1**2 / 2.d0
  m(4,5) = - hm2**2 / 2.d0
  m(4,6) = - hm3**2 / 2.d0
  m(4,7) = - hm4**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) = - hp1**3 / 6.d0
  m(5,3) =   0.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   hm2**3 / 6.d0
  m(5,6) =   hm3**3 / 6.d0
  m(5,7) =   hm4**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hp1**4 / 24.d0
  m(6,3) =   0.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) = - hm2**4 / 24.d0
  m(6,6) = - hm3**4 / 24.d0
  m(6,7) = - hm4**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) = - hp1**5 / 120.d0
  m(7,3) =   0.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   hm2**5 / 120.d0
  m(7,6) =   hm3**5 / 120.d0
  m(7,7) =   hm4**5 / 120.d0

 return

 end subroutine jmaxm1_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_jmaxm1_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_jmaxm1_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine jmaxm2_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hm3, hm2, hm1, dy

  dy = 1.d0
  ! 'normalized' length
  hp1 =       stf**3
  hp2 = hp1 + stf**4
  hm1 =       stf**2
  hm2 = hm1 + stf
  hm3 = hm2 + dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) = - hp2
  m(3,3) = - hp1
  m(3,4) =   0.d0
  m(3,5) =   hm1
  m(3,6) =   hm2
  m(3,7) =   hm3
 
  m(4,1) =   0.d0
  m(4,2) = - hp2**2 / 2.d0
  m(4,3) = - hp1**2 / 2.d0
  m(4,4) =   0.d0
  m(4,5) = - hm1**2 / 2.d0
  m(4,6) = - hm2**2 / 2.d0
  m(4,7) = - hm3**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) = - hp2**3 / 6.d0
  m(5,3) = - hp1**3 / 6.d0
  m(5,4) =   0.d0
  m(5,5) =   hm1**3 / 6.d0
  m(5,6) =   hm2**3 / 6.d0
  m(5,7) =   hm3**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hp2**4 / 24.d0
  m(6,3) = - hp1**4 / 24.d0
  m(6,4) =   0.d0
  m(6,5) = - hm1**4 / 24.d0
  m(6,6) = - hm2**4 / 24.d0
  m(6,7) = - hm3**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) = - hp2**5 / 120.d0
  m(7,3) = - hp1**5 / 120.d0
  m(7,4) =   0.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   hm2**5 / 120.d0
  m(7,7) =   hm3**5 / 120.d0

 return

 end subroutine jmaxm2_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_jmaxm2_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_jmaxm2_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine jmaxm3_point_wall_matrix_coef(m, n)
  
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hp3, hm2, hm1, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 =       stf**2
  hp2 = hp1 + stf**3
  hp3 = hp2 + stf**4
  hm1 = stf
  hm2 = hm1 + dy
  
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) = - hp3
  m(3,3) = - hp2
  m(3,4) = - hp1
  m(3,5) =   0.d0
  m(3,6) =   hm1
  m(3,7) =   hm2
 
  m(4,1) =   0.d0
  m(4,2) = - hp3**2 / 2.d0
  m(4,3) = - hp2**2 / 2.d0
  m(4,4) = - hp1**2 / 2.d0
  m(4,5) =   0.d0
  m(4,6) = - hm1**2 / 2.d0
  m(4,7) = - hm2**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) = - hp3**3 / 6.d0
  m(5,3) = - hp2**3 / 6.d0
  m(5,4) = - hp1**3 / 6.d0
  m(5,5) =   0.d0
  m(5,6) =   hm1**3 / 6.d0
  m(5,7) =   hm2**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hp3**4 / 24.d0
  m(6,3) = - hp2**4 / 24.d0
  m(6,4) = - hp1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hm1**4 / 24.d0
  m(6,7) = - hm2**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) = - hp3**5 / 120.d0
  m(7,3) = - hp2**5 / 120.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) =   hm1**5 / 120.d0
  m(7,7) =   hm2**5 / 120.d0

 return

 end subroutine jmaxm3_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_jmaxm3_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_jmaxm3_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine jmaxm4_point_wall_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hp1, hp2, hp3, hp4, hm1, dy
 
  dy = 1.d0
  ! 'normalized' length
  hp1 =       stf
  hp2 = hp1 + stf**2
  hp3 = hp2 + stf**3
  hp4 = hp3 + stf**4
  hm1 = dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) = - hp4
  m(3,3) = - hp3
  m(3,4) = - hp2
  m(3,5) = - hp1
  m(3,6) =   0.d0
  m(3,7) =   hm1
 
  m(4,1) =   0.d0
  m(4,2) = - hp4**2 / 2.d0
  m(4,3) = - hp3**2 / 2.d0
  m(4,4) = - hp2**2 / 2.d0
  m(4,5) = - hp1**2 / 2.d0
  m(4,6) =   0.d0
  m(4,7) = - hm1**2 / 2.d0
  
  m(5,1) =   0.d0
  m(5,2) = - hp4**3 / 6.d0
  m(5,3) = - hp3**3 / 6.d0
  m(5,4) = - hp2**3 / 6.d0
  m(5,5) = - hp1**3 / 6.d0
  m(5,6) =   0.d0
  m(5,7) =   hm1**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hp4**4 / 24.d0
  m(6,3) = - hp3**4 / 24.d0
  m(6,4) = - hp2**4 / 24.d0
  m(6,5) = - hp1**4 / 24.d0
  m(6,6) =   0.d0
  m(6,7) = - hm1**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) = - hp4**5 / 120.d0
  m(7,3) = - hp3**5 / 120.d0
  m(7,4) = - hp2**5 / 120.d0
  m(7,5) = - hp1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) =   hm1**5 / 120.d0

 return

 end subroutine jmaxm4_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine rhs_jmaxm4_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine rhs_jmaxm4_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine middle_first_point_wall_matrix_coef(m, n)
 
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
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) =   0.d0
  m(3,3) = - hp1
  m(3,4) = - hp2
  m(3,5) = - hp3
  m(3,6) = - hp4
  m(3,7) = - hp5
 
  m(4,1) =   0.d0
  m(4,2) =   0.d0
  m(4,3) = - hp1**2 / 2.d0
  m(4,4) = - hp2**2 / 2.d0
  m(4,5) = - hp3**2 / 2.d0
  m(4,6) = - hp4**2 / 2.d0
  m(4,7) = - hp5**2 / 2.d0

  m(5,1) =   0.d0
  m(5,2) =   0.d0
  m(5,3) = - hp1**3 / 6.d0
  m(5,4) = - hp2**3 / 6.d0
  m(5,5) = - hp3**3 / 6.d0
  m(5,6) = - hp4**3 / 6.d0
  m(5,7) = - hp5**3 / 6.d0

  m(6,1) =   0.d0
  m(6,2) =   0.d0
  m(6,3) = - hp1**4 / 24.d0
  m(6,4) = - hp2**4 / 24.d0
  m(6,5) = - hp3**4 / 24.d0
  m(6,6) = - hp4**4 / 24.d0
  m(6,7) = - hp5**4 / 24.d0

  m(7,1) =   0.d0
  m(7,2) =   0.d0
  m(7,3) = - hp1**5 / 120.d0
  m(7,4) = - hp2**5 / 120.d0
  m(7,5) = - hp3**5 / 120.d0
  m(7,6) = - hp4**5 / 120.d0
  m(7,7) = - hp5**5 / 120.d0

 return

 end subroutine middle_first_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_rhs_first_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine middle_rhs_first_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_second_point_wall_matrix_coef(m, n)
 
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
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0

  m(3,1) =   1.d0
  m(3,2) =   hm1
  m(3,3) =   0.d0
  m(3,4) = - hp1
  m(3,5) = - hp2
  m(3,6) = - hp3
  m(3,7) = - hp4
 
  m(4,1) =   0.d0
  m(4,2) = - hm1**2 / 2.d0
  m(4,3) =   0.d0
  m(4,4) = - hp1**2 / 2.d0
  m(4,5) = - hp2**2 / 2.d0
  m(4,6) = - hp3**2 / 2.d0
  m(4,7) = - hp4**2 / 2.d0

  m(5,1) =   0.d0
  m(5,2) =   hm1**3 / 6.d0
  m(5,3) =   0.d0
  m(5,4) = - hp1**3 / 6.d0
  m(5,5) = - hp2**3 / 6.d0
  m(5,6) = - hp3**3 / 6.d0
  m(5,7) = - hp4**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hm1**4 / 24.d0
  m(6,3) =   0.d0
  m(6,4) = - hp1**4 / 24.d0
  m(6,5) = - hp2**4 / 24.d0
  m(6,6) = - hp3**4 / 24.d0
  m(6,7) = - hp4**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hm1**5 / 120.d0
  m(7,3) =   0.d0
  m(7,4) = - hp1**5 / 120.d0
  m(7,5) = - hp2**5 / 120.d0
  m(7,6) = - hp3**5 / 120.d0
  m(7,7) = - hp4**5 / 120.d0

 return

 end subroutine middle_second_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_rhs_second_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine middle_rhs_second_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine middle_third_point_wall_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm2, hm1, hp1, hp2, hp3, dy

  dy = 1.d0
  ! 'normalized' length
  hm2 = dy  + stf    * dy
  hm1 =       stf    * dy
  hp1 =       stf**2 * dy
  hp2 = hp1 + stf**3 * dy
  hp3 = hp2 + stf**4 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) =   hm2
  m(3,3) =   hm1
  m(3,4) =   0.d0
  m(3,5) = - hp1
  m(3,6) = - hp2
  m(3,7) = - hp3
 
  m(4,1) =   0.d0
  m(4,2) = - hm2**2 / 2.d0
  m(4,3) = - hm1**2 / 2.d0
  m(4,4) =   0.d0
  m(4,5) = - hp1**2 / 2.d0
  m(4,6) = - hp2**2 / 2.d0
  m(4,7) = - hp3**2 / 2.d0

  m(5,1) =   0.d0
  m(5,2) =   hm2**3 / 6.d0
  m(5,3) =   hm1**3 / 6.d0
  m(5,4) =   0.d0
  m(5,5) = - hp1**3 / 6.d0
  m(5,6) = - hp2**3 / 6.d0
  m(5,7) = - hp3**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hm2**4 / 24.d0
  m(6,3) = - hm1**4 / 24.d0
  m(6,4) =   0.d0
  m(6,5) = - hp1**4 / 24.d0
  m(6,6) = - hp2**4 / 24.d0
  m(6,7) = - hp3**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hm2**5 / 120.d0
  m(7,3) =   hm1**5 / 120.d0
  m(7,4) =   0.d0
  m(7,5) = - hp1**5 / 120.d0
  m(7,6) = - hp2**5 / 120.d0
  m(7,7) = - hp3**5 / 120.d0

 return

 end subroutine middle_third_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_rhs_third_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine middle_rhs_third_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine middle_fourth_point_wall_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm3, hm2, hm1, hp1, hp2, dy
 
  dy = 1.d0
  ! 'normalized' length
  hm1 =       stf**2 * dy
  hm2 = hm1 + stf    * dy
  hm3 = hm2 + dy
  hp1 =       stf**3 * dy
  hp2 = hp1 + stf**4 * dy
 
  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0
 
  m(3,1) =   1.d0
  m(3,2) =   hm3
  m(3,3) =   hm2
  m(3,4) =   hm1
  m(3,5) =   0.d0
  m(3,6) = - hp1
  m(3,7) = - hp2
 
  m(4,1) =   0.d0
  m(4,2) = - hm3**2 / 2.d0
  m(4,3) = - hm2**2 / 2.d0
  m(4,4) = - hm1**2 / 2.d0
  m(4,5) =   0.d0
  m(4,6) = - hp1**2 / 2.d0
  m(4,7) = - hp2**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) =   hm3**3 / 6.d0
  m(5,3) =   hm2**3 / 6.d0
  m(5,4) =   hm1**3 / 6.d0
  m(5,5) =   0.d0
  m(5,6) = - hp1**3 / 6.d0
  m(5,7) = - hp2**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hm3**4 / 24.d0
  m(6,3) = - hm2**4 / 24.d0
  m(6,4) = - hm1**4 / 24.d0
  m(6,5) =   0.d0
  m(6,6) = - hp1**4 / 24.d0
  m(6,7) = - hp2**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hm3**5 / 120.d0
  m(7,3) =   hm2**5 / 120.d0
  m(7,4) =   hm1**5 / 120.d0
  m(7,5) =   0.d0
  m(7,6) = - hp1**5 / 120.d0
  m(7,7) = - hp2**5 / 120.d0

 return

 end subroutine middle_fourth_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_rhs_fourth_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine middle_rhs_fourth_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_fifth_point_wall_matrix_coef(m, n)
 
  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(out) :: m(n,n)
  real(kind=8)              :: hm4, hm3, hm2, hm1, hp1, dy
 
  dy = 1.d0
  ! 'normalized' length
  hm1 =       stf**3 * dy
  hm2 = hm1 + stf**2 * dy
  hm3 = hm2 + stf    * dy
  hm4 = hm3 + dy
  hp1 =       stf**4 * dy

  ! matrix definition
  m(1,1) =   1.d0
  m(1,2) =   0.d0
  m(1,3) =   0.d0
  m(1,4) =   0.d0
  m(1,5) =   0.d0
  m(1,6) =   0.d0
  m(1,7) =   0.d0
 
  m(2,1) =   0.d0
  m(2,2) = - 1.d0
  m(2,3) = - 1.d0
  m(2,4) = - 1.d0
  m(2,5) = - 1.d0
  m(2,6) = - 1.d0
  m(2,7) = - 1.d0

  m(3,1) =   1.d0
  m(3,2) =   hm4
  m(3,3) =   hm3
  m(3,4) =   hm2
  m(3,5) =   hm1
  m(3,6) =   0.d0
  m(3,7) = - hp1

  m(4,1) =   0.d0
  m(4,2) = - hm4**2 / 2.d0
  m(4,3) = - hm3**2 / 2.d0
  m(4,4) = - hm2**2 / 2.d0
  m(4,5) = - hm1**2 / 2.d0
  m(4,6) =   0.d0
  m(4,7) = - hp1**2 / 2.d0
 
  m(5,1) =   0.d0
  m(5,2) =   hm4**3 / 6.d0
  m(5,3) =   hm3**3 / 6.d0
  m(5,4) =   hm2**3 / 6.d0
  m(5,5) =   hm1**3 / 6.d0
  m(5,6) =   0.d0
  m(5,7) = - hp1**3 / 6.d0
 
  m(6,1) =   0.d0
  m(6,2) = - hm4**4 / 24.d0
  m(6,3) = - hm3**4 / 24.d0
  m(6,4) = - hm2**4 / 24.d0
  m(6,5) = - hm1**4 / 24.d0
  m(6,6) =   0.d0
  m(6,7) = - hp1**4 / 24.d0
 
  m(7,1) =   0.d0
  m(7,2) =   hm4**5 / 120.d0
  m(7,3) =   hm3**5 / 120.d0
  m(7,4) =   hm2**5 / 120.d0
  m(7,5) =   hm1**5 / 120.d0
  m(7,6) =   0.d0
  m(7,7) = - hp1**5 / 120.d0

 return

 end subroutine middle_fifth_point_wall_matrix_coef
 
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine middle_rhs_fifth_point_wall_matrix_coef(b, n)

  ! right-hand side of first derivative coef.
  implicit none
  integer, intent(in)       :: n
  integer                   :: i
  real(kind=8), intent(out) :: b(n)
 
  do i = 1, n
    b(i) = 0.d0
  end do
  b(1) = 1.d0
 
  return

 end subroutine middle_rhs_fifth_point_wall_matrix_coef

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine at_wall_coef_generation
 
  implicit none
  real(kind=8) :: fp_w_c(7,7), fp_w_b(7), fp_w_coef(7)
  real(kind=8) :: sp_w_c(7,7), sp_w_b(7), sp_w_coef(7)
  real(kind=8) :: tp_w_c(7,7), tp_w_b(7), tp_w_coef(7)
  real(kind=8) :: fourp_w_c(7,7), fourp_w_b(7), fourp_w_coef(7)
  real(kind=8) :: fivep_w_c(7,7), fivep_w_b(7), fivep_w_coef(7)
 
  real(kind=8) :: jmaxp_w_c(7,7), jmaxp_w_b(7), jmaxp_w_coef(7)
  real(kind=8) :: jmaxm1p_w_c(7,7), jmaxm1p_w_b(7), jmaxm1p_w_coef(7)
  real(kind=8) :: jmaxm2p_w_c(7,7), jmaxm2p_w_b(7), jmaxm2p_w_coef(7)
  real(kind=8) :: jmaxm3p_w_c(7,7), jmaxm3p_w_b(7), jmaxm3p_w_coef(7)
  real(kind=8) :: jmaxm4p_w_c(7,7), jmaxm4p_w_b(7), jmaxm4p_w_coef(7)
 
  real(kind=8) :: m_fp_w_c(7,7), m_fp_w_b(7), m_fp_w_coef(7)
  real(kind=8) :: m_sp_w_c(7,7), m_sp_w_b(7), m_sp_w_coef(7)
  real(kind=8) :: m_tp_w_c(7,7), m_tp_w_b(7), m_tp_w_coef(7)
  real(kind=8) :: m_fourp_w_c(7,7), m_fourp_w_b(7), m_fourp_w_coef(7)
  real(kind=8) :: m_fivep_w_c(7,7), m_fivep_w_b(7), m_fivep_w_coef(7)

  real(kind=8) :: m(5,7)
  real(kind=8) :: mm(5,5),  b(5),  x(5)
  real(kind=8) :: mm2(4,4), b2(4), x2(4)
 
  real(kind=8) :: mf(5,7)
  real(kind=8) :: mmf(5,5), bf(5), xf(5)
 
  real(kind=8) :: xm(5)
 
! real(kind=8) :: sp_integ_coef(7)
! real(kind=8) :: mp_integ_coef(8)
! real(kind=8) :: pp_integ_coef(7)
! real(kind=8) :: lp_integ_coef(7)

  ! First point coefficients
  call first_point_wall_matrix_coef(fp_w_c, 7)
  call rhs_first_point_wall_matrix_coef(fp_w_b, 7)
  call ludecomp(fp_w_c, fp_w_b, fp_w_coef, 7)
 
  ! Second point coefficients
  call second_point_wall_matrix_coef(sp_w_c, 7)
  call rhs_second_point_wall_matrix_coef(sp_w_b, 7)
  call ludecomp(sp_w_c, sp_w_b, sp_w_coef, 7)
 
 ! Third point coefficients
  call third_point_wall_matrix_coef(tp_w_c, 7)
  call rhs_third_point_wall_matrix_coef(tp_w_b, 7)
  call ludecomp(tp_w_c, tp_w_b, tp_w_coef, 7)
 
 ! Fourth point coefficients
  call fourth_point_wall_matrix_coef(fourp_w_c, 7)
  call rhs_fourth_point_wall_matrix_coef(fourp_w_b, 7)
  call ludecomp(fourp_w_c, fourp_w_b, fourp_w_coef, 7)

 ! Fifth point coefficients
  call fifth_point_wall_matrix_coef(fivep_w_c, 7)
  call rhs_fifth_point_wall_matrix_coef(fivep_w_b, 7)
  call ludecomp(fivep_w_c, fivep_w_b, fivep_w_coef, 7)

 ! Last point coefficients
  call jmax_point_wall_matrix_coef(jmaxp_w_c, 7)
  call rhs_jmax_point_wall_matrix_coef(jmaxp_w_b, 7)
  call ludecomp(jmaxp_w_c, jmaxp_w_b, jmaxp_w_coef, 7)
 
  call jmaxm1_point_wall_matrix_coef(jmaxm1p_w_c, 7)
  call rhs_jmaxm1_point_wall_matrix_coef(jmaxm1p_w_b, 7)
  call ludecomp(jmaxm1p_w_c, jmaxm1p_w_b, jmaxm1p_w_coef, 7)
 
  call jmaxm2_point_wall_matrix_coef(jmaxm2p_w_c, 7)
  call rhs_jmaxm2_point_wall_matrix_coef(jmaxm2p_w_b, 7)
  call ludecomp(jmaxm2p_w_c, jmaxm2p_w_b, jmaxm2p_w_coef, 7)
 
  call jmaxm3_point_wall_matrix_coef(jmaxm3p_w_c, 7)
  call rhs_jmaxm3_point_wall_matrix_coef(jmaxm3p_w_b, 7)
  call ludecomp(jmaxm3p_w_c, jmaxm3p_w_b, jmaxm3p_w_coef, 7)
 
  call jmaxm4_point_wall_matrix_coef(jmaxm4p_w_c, 7)
  call rhs_jmaxm4_point_wall_matrix_coef(jmaxm4p_w_b, 7)
  call ludecomp(jmaxm4p_w_c, jmaxm4p_w_b, jmaxm4p_w_coef, 7)
 
  ! Middle first point coefficients
  call middle_first_point_wall_matrix_coef(m_fp_w_c, 7)
  call middle_rhs_first_point_wall_matrix_coef(m_fp_w_b, 7)
  call ludecomp(m_fp_w_c, m_fp_w_b, m_fp_w_coef, 7)
 
  ! Middle second point coefficients
  call middle_second_point_wall_matrix_coef(m_sp_w_c, 7)
  call middle_rhs_second_point_wall_matrix_coef(m_sp_w_b, 7)
  call ludecomp(m_sp_w_c, m_sp_w_b, m_sp_w_coef, 7)
 
  ! Middle third point coefficients
  call middle_third_point_wall_matrix_coef(m_tp_w_c, 7)
  call middle_rhs_third_point_wall_matrix_coef(m_tp_w_b, 7)
  call ludecomp(m_tp_w_c, m_tp_w_b, m_tp_w_coef, 7)
 
  ! Middle fourth point coefficients
  call middle_fourth_point_wall_matrix_coef(m_fourp_w_c, 7)
  call middle_rhs_fourth_point_wall_matrix_coef(m_fourp_w_b, 7)
  call ludecomp(m_fourp_w_c, m_fourp_w_b, m_fourp_w_coef, 7)

  ! Middle fifth point coefficients
  call middle_fifth_point_wall_matrix_coef(m_fivep_w_c, 7)
  call middle_rhs_fifth_point_wall_matrix_coef(m_fivep_w_b, 7)
  call ludecomp(m_fivep_w_c, m_fivep_w_b, m_fivep_w_coef, 7)
 
  m(1,:) = fp_w_coef(:)
  m(2,:) = sp_w_coef(:)
  m(3,:) = tp_w_coef(:)
  m(4,:) = fourp_w_coef(:)
  m(5,:) = fivep_w_coef(:)

  mm(1,:) = 0.d0
  mm(1,1) = 1.d0
  mm(2,:) = m(:,4)
  mm(3,:) = m(:,5)
  mm(4,:) = m(:,6)
  mm(5,:) = m(:,7)

  b(:) = 0.d0
  b(1) = 251.d0
 
  call ludecomp(mm, b, x, 5)
  sp_integ_coef(1)   = - dot_product(m(:,2),x)
  sp_integ_coef(2)   = - dot_product(m(:,3),x)
  sp_integ_coef(3:7) = x
! print *, 'Second point coef:'
! print *, sp_integ_coef

  ! Middle app calculation
  m(1,:) = m_fp_w_coef(:)
  m(2,:) = m_sp_w_coef(:)
  m(3,:) = m_tp_w_coef(:)
  m(4,:) = m_fourp_w_coef(:)
  m(5,:) = m_fivep_w_coef(:)

  mm(1,:) = 0.d0
  mm(1,1) = 1.d0
  mm(2,:) = 0.d0
  mm(2,2) = 1.d0
  mm(3,:) = m(:,5)
  mm(4,:) = m(:,6)
  mm(5,:) = m(:,7)

  b(:) = 0.d0
  b(1) = 281.d0
  b(2) = 2056.d0
 
  call ludecomp(mm, b, xm, 5)

  mp_integ_coef(1)   = - dot_product(m(:,2),xm)
  mp_integ_coef(2)   = - dot_product(m(:,3),xm)
  mp_integ_coef(3)   = - dot_product(m(:,4),xm)
  mp_integ_coef(4:8) = xm
! print *, 'Middle point coef:'
! print *, mp_integ_coef

  ! Penultimate point calculation
  mm2(1,:) = 0.d0
  mm2(1,1) = 1.d0
  mm2(2,:) = m(1:4,5)
  mm2(3,:) = m(1:4,6)
  mm2(4,:) = m(1:4,7)
 
  b2(:) = 0.d0
  b2(1) = 10.d0
 
  call ludecomp(mm2, b2, x2, 4)
  pp_integ_coef(1)   = - dot_product(m(1:4,2),x2)
  pp_integ_coef(2)   = - dot_product(m(1:4,3),x2)
  pp_integ_coef(3)   = - dot_product(m(1:4,4),x2)
  pp_integ_coef(4:7) = x2
! print *, 'penultimate point coef:'
! print *, pp_integ_coef

  ! Last point calculation
  mf(1,:) = jmaxp_w_coef(:)
  mf(2,:) = jmaxm1p_w_coef(:)
  mf(3,:) = jmaxm2p_w_coef(:)
  mf(4,:) = jmaxm3p_w_coef(:)
  mf(5,:) = jmaxm4p_w_coef(:)
 
  mmf(1,:) = 0.d0
  mmf(1,1) = 1.d0
  mmf(2,:) = mf(:,4)
  mmf(3,:) = mf(:,5)
  mmf(4,:) = mf(:,6)
  mmf(5,:) = mf(:,7)

  bf(:) = 0.d0
  bf(1) = 251.d0
 
  call ludecomp(mmf, bf, xf, 5)
 
  lp_integ_coef(1)   = - dot_product(mf(:,2),xf)
  lp_integ_coef(2)   = - dot_product(mf(:,3),xf)
  lp_integ_coef(3:7) = xf
! print *, 'Last point coef:'
! print *, lp_integ_coef
! 
! print *, xf

 end subroutine at_wall_coef_generation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

end module atwall_calculation

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
