cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                       begin of pentadiagonal solvers                         c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine bandy5(a, al, indx)

      ! solve the LHS of one pentadiagonal matrix in y direction
      ! the a is the input and a, al and indx are outputs
      implicit none
      include '../par.for'
      integer m1, indx(jmax), i, j, k, l, mm
      real*8 d, a(jmax,5), al(jmax,5), tiny, dum
      parameter (tiny = 1.d-20)

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
      do k = 1, jmax
        dum = a(k,1)
        i   = k
        if (l.lt.jmax) l = l + 1
        indx(k) = i
        if (dum .eq. 0.d0) a(k,1) = tiny
        if (i .ne. k) then
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
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine banbky5(a, al, indx, rhs)

      ! solve the the pentadiagonal matrix in y direction
      ! the terms a and al come from the subroutine bandy5
      ! the rhs variable is the input and at the end
      ! is the result of the solved problem
      implicit none
      include '../par.for'
      integer indx(jmax), i, k, l
      real*8 a(jmax,5), al(jmax,5), rhs(jmax), dum

      l = 2
      do k = 1, jmax
        i = indx(k)
        if (i .ne. k) then
          dum    = rhs(k)
          rhs(k) = rhs(i)
          rhs(i) = dum
        end if
        if (l .lt. jmax) l = l + 1
        do i = k + 1, l
          rhs(i) = rhs(i) - al(k,i-k) * rhs(k)
        end do
      end do
      l = 1
      do i = jmax, 1, -1
        dum = rhs(i)
        do k = 2, l
          dum = dum - a(i,k) * rhs(k+i-1)
        end do
        rhs(i) = dum / a(i,1)
        if (l .lt. 5) l = l + 1
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine band5_poi(a, n, al, indy)

      ! solve the LHS of the pentadiagonal matrix
      ! n is the size of the matrix
      ! used for uy poisson solver subroutine
      implicit none
      include '../par.for'
      integer m1, n, indy(jmax), i, j, k, l, mm
      real*8 d, a(jmax,5), al(jmax,5), TINY, dum
      parameter (TINY = 1.e-20)

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
        if (l .lt. n) l = l + 1
        indy(k) = i
        if (dum .eq. 0.d0) a(k,1) = TINY
        if (i .ne. k) then
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
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine banbk5_poi(a, n, al, indy, rhs)

      ! solve the the pentadiagonal matrix the term a and al comes from
      ! subroutine band5_poi, the rhs variable is the input and at the
      ! end is the result of matrix
      ! n is the size of the matrix
      ! used for uy poisson solver subroutine
      implicit none
      include '../par.for'
      integer n, i, k, l, indy(jmax)
      real*8 a(jmax,5), al(jmax,5), rhs(jmax), dum

      l = 2
      do k = 1, n
        i = indy(k)
        if (i .ne. k) then
          dum    = rhs(k)
          rhs(k) = rhs(i)
          rhs(i) = dum
        end if
        if (l .lt. n) l = l + 1
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
        if (l .lt. 5) l = l + 1
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                        end of pentadiagonal solvers                          c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
