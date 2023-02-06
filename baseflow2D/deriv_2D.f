cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         derivative calculations                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine der(df, fc, delta, n)

      ! first derivatives calculation in x direction
      implicit none
      integer n
      real*8 a1x(600), b1x(600), c1x(600), fc(600), df(600), delta
      common/der1x/ a1x, b1x, c1x

      call coef2(a1x, b1x, c1x, n)
      call rhs(fc, df, delta, n)
      call trid(a1x, b1x, c1x, df, n)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coef2(a, b, c, n)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      integer i, n
      real*8 a(600), b(600), c(600)

      a(1)   = 0.d0
      b(1)   = 1.d0
      c(1)   = 0.d0

      a(2)   = 1.d0
      b(2)   = 6.d0
      c(2)   = 2.d0

      do i = 3, n-2
        a(i) = 1.d0
        b(i) = 3.d0
        c(i) = 1.d0
      end do

      a(n-1) = 2.d0
      b(n-1) = 6.d0
      c(n-1) = 1.d0

      a(n)   = 4.d0
      b(n)   = 1.d0
      c(n)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhs(fc, df, delta, n)

      ! RHS for the first derivative calculation in x direction
      implicit none
      integer i, n
      real*8 df(600), fc(600), delta

      df(1) = ( - 25.d0 * fc(1) + 48.d0 * fc(2)
     &          - 36.d0 * fc(3) + 16.d0 * fc(4)
     &          -  3.d0 * fc(5) ) / ( 12.d0 * delta )

      df(2) = ( - 406.d0 * fc(1) - 300.d0 * fc(2)
     &          + 760.d0 * fc(3) -  80.d0 * fc(4)
     &          +  30.d0 * fc(5) -   4.d0 * fc(6) )
     &      /   ( 120.d0 * delta )

      do i = 3, n-2
        df(i) = (           fc(i+2) - fc(i-2)
     &          + 28.d0 * ( fc(i+1) - fc(i-1) ) )
     &        / ( 12.d0 * delta )
      end do

      df(n-1) = ( - 406.d0 * fc(n)
     &            - 300.d0 * fc(n-1)
     &            + 760.d0 * fc(n-2)
     &            -  80.d0 * fc(n-3)
     &            +  30.d0 * fc(n-4)
     &            -   4.d0 * fc(n-5) )
     &        / ( - 120.d0 * delta )

      df(n) = ( - 74.d0 * fc(n)
     &          + 16.d0 * fc(n-1)
     &          + 72.d0 * fc(n-2)
     &          - 16.d0 * fc(n-3)
     &          +  2.d0 * fc(n-4) )
     &      / ( - 24.d0 * delta )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine trid(a, b, c, rhs, n)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      integer i, n
      real*8 a(600), b(600), c(600), gam(600), bet, rhs(600), u(600)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do i = 2, n
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
        u(i)   = ( rhs(i) - a(i) * u(i-1) ) / bet
      end do
      do i = n - 1, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do
      do i = 1, n
        rhs(i) = u(i)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                      end of derivative calculations                          c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
