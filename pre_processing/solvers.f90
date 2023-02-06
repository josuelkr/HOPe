!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

module solvers

 implicit none

 contains

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine tridiagonal(a, b, c, r, n)

  implicit none
  integer, intent(in)            :: n
  real(kind=8), intent(in)       :: a(n), b(n), c(n)
  real(kind=8), intent(inout)    :: r(n)
 
  integer                        :: j
  real(kind=8)                   :: gam(n), u(n), bet
 
  bet  = b(1)
  u(1) = r(1) / bet
  do j = 2, n
   gam(j) = c(j-1) / bet
   bet    = b(j) - a(j) * gam(j)
   u(j)   = ( r(j) - a(j) * u(j-1) ) / bet
  end do
  do j = n - 1, 1, -1
   u(j) = u(j) - gam(j+1) * u(j+1)
  end do
 
  r = u ! Vector operation
 
  return 

 end subroutine tridiagonal

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine tridiagonal_hp(a, b, c, r, n)

  implicit none
  integer, intent(in)            :: n
  real(kind=16), intent(in)      :: a(n), b(n), c(n)
  real(kind=16), intent(inout)   :: r(n)
 
  integer                        :: j
  real(kind=16)                  :: gam(n), u(n), bet
 
  bet  = b(1)
  u(1) = r(1) / bet
  do j = 2, n
   gam(j) = c(j-1) / bet
   bet    = b(j) - a(j) * gam(j)
   u(j)   = ( r(j) - a(j) * u(j-1) ) / bet
  end do
  do j = n - 1, 1, -1
   u(j) = u(j) - gam(j+1) * u(j+1)
  end do
 
  r = u ! Vector operation
 
  return 

 end subroutine tridiagonal_hp

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

 subroutine ludecomp(a, b, x, n)

  ! LU solver
 
  implicit none
  integer, intent(in)           :: n
  real(kind=8), intent(in)      :: a(n,n)
  real(kind=8), intent(in)      :: b(n)
  real(kind=8), intent(out)     :: x(n)
  real(kind=8)                  :: l(n,n),u(n,n),y(n),soma
  integer                       :: i,j,it,k
 
  ! LU decomposition
  do it = 1, n
    i = it
    do j = it, n
      soma = 0.d0
      do k = 1, i - 1
       soma = soma + l(i,k) * u(k,j)
      end do
      u(i,j) = a(i,j) - soma
    end do
    j = it
    do i = it, n
      soma = 0.d0
      do k = 1, j - 1
       soma = soma + l(i,k) * u(k,j)
      end do
      l(i,j) = (a(i,j) - soma) / u(j,j)
    end do
  end do
  do i = 1, n
    l(i,i) = 1.d0
  end do
 
  ! lower system resolution
  y(1) = b(1) / l(1,1)
  do i = 2, n
    soma = 0.d0
    do j = 1, i - 1
      soma = soma + l(i,j) * y(j)
    end do
    y(i) = (b(i) - soma) / l(i,i)
  end do
 
  ! upper system resolution
  x(n) = y(n) / u(n,n)
  do i = n-1, 1, -1
    soma = 0.d0
    do j = i + 1, n
      soma = soma + u(i,j) * x(j)
    end do
    x(i) = (y(i) - soma) / u(i,i)
    if (dabs(x(i)) < 1d-14) x(i) = 0.d0
  end do
 
   return

  end subroutine ludecomp

!---------------------------------------------------------------------------------
!---------------------------------------------------------------------------------

 subroutine ludecomp_hp(a, b, x, n)

  ! LU solver
 
  implicit none
  integer, intent(in)           :: n
  real(kind=16), intent(in)     :: a(n,n)
  real(kind=16), intent(in)     :: b(n)
  real(kind=16), intent(out)    :: x(n)
  real(kind=16)                 :: l(n,n),u(n,n),y(n),soma
  integer                       :: i,j,it,k

  ! LU decomposition
  do it = 1, n
    i = it
    do j = it, n
      soma = 0.d0
      do k = 1, i - 1
       soma = soma + l(i,k) * u(k,j)
      end do
      u(i,j) = a(i,j) - soma
    end do
    j = it
    do i = it, n
      soma = 0.d0
      do k = 1, j - 1
       soma = soma + l(i,k) * u(k,j)
      end do
      l(i,j) = (a(i,j) - soma) / u(j,j)
    end do
  end do
  do i = 1, n
    l(i,i) = 1.d0
  end do
 
  ! lower system resolution
  y(1) = b(1) / l(1,1)
  do i = 2, n
    soma = 0.d0
    do j = 1, i - 1
      soma = soma + l(i,j) * y(j)
    end do
    y(i) = (b(i) - soma) / l(i,i)
  end do
 
  ! upper system resolution
  x(n) = y(n) / u(n,n)
  do i = n-1, 1, -1
    soma = 0.d0
    do j = i + 1, n
      soma = soma + u(i,j) * x(j)
    end do
    x(i) = (y(i) - soma) / u(i,i)
    if (abs(x(i)) < 1d-14) x(i) = 0.d0
  end do
 
   return

  end subroutine ludecomp_hp

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine bandy5(a, al, indx, n)

  ! solve the LHS of one pentadiagonal matrix in y direction
  ! the a and n variables are inputs and a, al and indx are output
 
  implicit none
  integer, intent(in)           :: n
  real(kind=8), intent(inout)   :: a(n,5)
  real(kind=8), intent(out)     :: al(n,5)
  integer, intent(out)          :: indx(n)
 
  integer                       :: m1, i, j, k, l, mm
  real(kind=8)                  :: d, dum
  real(kind=8), parameter       :: tiny=1.d-20
 
   m1 = 2
   mm = 5
   l  = m1
   do i = 1, m1
     do j = m1 + 2 - i, mm
       a(i,j-l) = a(i,j)
     end do
     l = l - 1
     do j = mm - l, mm
       a(i,j) = 0.d0
     end do
   end do
   d = 1.d0
   l = m1
   do k = 1, n
     dum = a(k,1)
     i   = k
     if (l < n) l = l + 1
     do j = k + 1, l
       if (dabs(a(j,1)) > dabs(dum)) then
         dum = a(j,1)
         i   = j
       end if
     end do
     indx(k) = i
     if (dum == 0.d0) a(k,1) = tiny
     if (i /= k) then
       d = - d
       do j = 1, mm
         dum    = a(k,j)
         a(k,j) = a(i,j)
         a(i,j) = dum
       end do
     end if
     do i = k + 1, l
       dum       = a(i,1) / a(k,1)
       al(k,i-k) = dum
       do j = 2, mm
         a(i,j-1) = a(i,j) - dum * a(k,j)
       end do
       a(i,mm) = 0.d0
     end do
   end do
 
  return

 end subroutine bandy5

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine banbky5(a, al, indx, rhs, n)

  ! solve the the pentadiagonal matrix in y direction
  ! the terms a and al come from the subroutine bandy5
  ! the rhs variable is the input and at the end
  ! is the result of the solved problem
 
  implicit none
  integer, intent(in)            :: n
  integer, intent(in)            :: indx(n)
  real(kind=8), intent(in)       :: a(n,5), al(n,5)
  real(kind=8), intent(inout)    :: rhs(n)
 
  integer                        :: i, k, l
  real(kind=8)                   :: dum
 
  l = 2
  do k = 1, n
    i = indx(k)
    if (i /= k) then
      dum    = rhs(k)
      rhs(k) = rhs(i)
      rhs(i) = dum
    end if
    if (l < n) l = l + 1
    do i = k + 1, l
      rhs(i) = rhs(i) - al(k,i-k) * rhs(k)
    end do
  end do
  l = 1
  do i = n, 1, -1
    dum = rhs(i)
    do k = 2, l
      dum = dum - a(i,k) * rhs(k+i-1)
    end do
    rhs(i) = dum / a(i,1)
    if (l < 5) l = l + 1
  end do
 
  return

 end subroutine banbky5

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
 
 subroutine alternative_ludecomp(a, lu, n)

  implicit none
  integer, intent(in)       :: n
  real(kind=8), intent(in)  :: a(n,5)
  real(kind=8), intent(out) :: lu(n,5)
  integer                   :: i, j, k, it, jj, ii, kk
  real(kind=8)              :: soma
 
  do i = 1, n - 2
    jj = i - 1
    if (jj > 2) jj = 2
    do j = 3, 5
      kk   = 2
      soma = 0.d0
      do k = 1, jj
        soma = soma + lu(i,kk) * lu(i-k,j+k)
        kk   = kk - 1
      end do
      jj = jj - 1
      lu(i,j) = a(i,j) - soma
    end do
    j  = 2
    ii = i - 1
    if (ii > 1) ii = 1
    do it = i + 1, i + 2
      soma = 0.d0
      do k = 1, ii
        soma = lu(it,k) * lu(i-1,j+2)
      end do
      lu(it,j) = ( a(it,j) - soma ) / lu(i,3)
      ii = ii - 1
      j  = j - 1
    end do
  end do
  i = n - 1
  lu(i,3)   =   a(i,3)   - lu(i,2)   * lu(i-1,4)   - lu(i,1) * lu(i-2,5) 
  lu(i,4)   =   a(i,4)   - lu(i,2)   * lu(i-1,5)
  lu(i+1,2) = ( a(i+1,2) - lu(i+1,1) * lu(i-1,4) ) / lu(i,3)
  i = n
  lu(i,3)   =   a(i,3)   - lu(i,2)   * lu(i-1,4)   - lu(i,1) * lu(i-2,5) 
 
  return

 end subroutine alternative_ludecomp

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

 subroutine alternative_lusolver(lu, rhs, n)

  implicit none
  integer, intent(in)         :: n
  real(kind=8), intent(in)    :: lu(n,5)
  real(kind=8), intent(inout) :: rhs(n)
  integer                     :: j
  real(kind=8)                :: raux(3), y(n)
 
!  begin of solver Ly = rhs
   y(1) = rhs(1)
   y(2) = rhs(2) - lu(2,2) * y(1)
   do j = 3 , n
    y(j) = rhs(j) - lu(j,1) * y(j-2) - lu(j,2) * y(j-1)
   end do
!  end of solver Ly = rhs

!  begin of solver U(rhs) = y
   rhs(n)   = y(n) / lu(n,3)
   rhs(n-1) = (y(n-1) - lu(n-1,4) * rhs(n)) / lu(n-1,3)
   do j = n - 2, 1, -1
    rhs(j) = ( y(j) - lu(j,4) * rhs(j+1) - lu(j,5) * rhs(j+2) ) / lu(j,3)
   end do
!  end of solver U(rhs) = y

  return

 end subroutine alternative_lusolver

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

end module solvers

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
