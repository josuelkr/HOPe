cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                           derivative calculations                            c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivs_k

      ! mounts the tri-diagonal LHS for derivative calculations
      implicit none
      include '../par.for'
      real*8  a1x(ptsx),  b1x(ptsx),  c1x(ptsx)
      real*8  a2x(ptsx),  b2x(ptsx),  c2x(ptsx)
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      real*8  a2y(jmax),  b2y(jmax),  c2y(jmax)
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      common/der1x/  a1x,  b1x,  c1x
      common/der2x/  a2x,  b2x,  c2x
      common/der1y/  a1y,  b1y,  c1y
      common/der2y/  a2y,  b2y,  c2y
      common/der1fv/ a1fv, b1fv, c1fv

      call coefx(a1x, b1x, c1x)
      call coeffx(a2x, b2x, c2x)
      call coef(a1y, b1y, c1y)
      call coeffv(a1fv, b1fv, c1fv)
      call coeff(a2y, b2y, c2y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxue(ddx, fc)

      ! first derivatives calculation in x direction (parallel)
      implicit none
      include '../par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx)
      real*8 fc(ptsx), ddx(ptsx)
      common/der1x/ a1x, b1x, c1x

      call rhsxue(fc, ddx)
      call tridparxue(a1x, b1x, c1x, ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparx(ddx, fc)

      ! first derivatives calculation in x direction (parallel)
      implicit none
      include '../par.for'
      real*8 a1x(ptsx), b1x(ptsx), c1x(ptsx)
      real*8 fc(ptsx,jmax), ddx(ptsx,jmax)
      common/der1x/ a1x, b1x, c1x

      call rhsx(fc, ddx)
      call tridparx(a1x, b1x, c1x, ddx)
c     call tridseqx(a1x, b1x, c1x, ddx)
c     call boundary_exchange_derivs(ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derparxx(d2dx2, fc)

      ! second derivative calculations in x direction (parallel)
      implicit none
      include '../par.for'
      real*8 a2x(ptsx), b2x(ptsx), c2x(ptsx)
      real*8 fc(ptsx,jmax), d2dx2(ptsx,jmax)
      common/der2x/ a2x, b2x, c2x

      call rhsxx(fc, d2dx2)
      call tridparx(a2x, b2x, c2x, d2dx2)
c     call tridseqx(a2x, b2x, c2x, d2dx2)
c     call boundary_exchange_derivs(d2dx2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dery(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include '../par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      real*8 fc(ptsx,jmax), ddy(ptsx,jmax)
      common/der1y/ a1y, b1y, c1y

      call rhsy(fc, ddy)
      call tridy(a1y, b1y, c1y, ddy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryfv(ddy, fc)

      ! first derivative calculation in y direction for uy
      implicit none
      include '../par.for'
      real*8 a1fv(jmax), b1fv(jmax), c1fv(jmax)
      real*8 fc(ptsx,jmax), ddy(ptsx,jmax)
      common/der1fv/ a1fv, b1fv, c1fv

      call rhsyfv(fc, ddy)
      call tridy(a1fv, b1fv, c1fv, ddy)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryy(d2dy2, fc)

      ! second derivative calculation in y direction
      implicit none
      include '../par.for'
      real*8 a2y(jmax), b2y(jmax), c2y(jmax)
      real*8 fc(ptsx,jmax), d2dy2(ptsx,jmax)
      common/der2y/ a2y, b2y, c2y

      call rhsyy(fc, d2dy2)
      call tridy(a2y, b2y, c2y, d2dy2)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxue(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include '../par.for'
      integer i
      real*8 rhs(ptsx), fc(ptsx)

      rhs(1) = ( - 25.d0 * fc(1) + 48.d0 * fc(2)
     &           - 36.d0 * fc(3) + 16.d0 * fc(4)
     &           -  3.d0 * fc(5) ) / ( 12.d0 * dx )
      
      rhs(2) = ( - 406.d0 * fc(1) - 300.d0 * fc(2)
     &           + 760.d0 * fc(3) -  80.d0 * fc(4)
     &           +  30.d0 * fc(5) -   4.d0 * fc(6) )
     &       /   ( 120.d0 * dx )
      
      do i = 3, ptsx-2
        rhs(i) = (           fc(i+2) - fc(i-2)
     &           + 28.d0 * ( fc(i+1) - fc(i-1) ) )
     &         / ( 12.d0 * dx )
      end do
      
      rhs(ptsx-1) = ( - 406.d0 * fc(ptsx)
     &                - 300.d0 * fc(ptsx-1)
     &                + 760.d0 * fc(ptsx-2)
     &                -  80.d0 * fc(ptsx-3)
     &                +  30.d0 * fc(ptsx-4)
     &                -   4.d0 * fc(ptsx-5) )
     &            / ( - 120.d0 * dx )
      
      rhs(ptsx) = ( - 74.d0 * fc(ptsx)
     &              + 16.d0 * fc(ptsx-1)
     &              + 72.d0 * fc(ptsx-2)
     &              - 16.d0 * fc(ptsx-3)
     &              +  2.d0 * fc(ptsx-4) )
     &          / ( - 24.d0 * dx )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsx(fc, rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include '../par.for'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do j = 1, jmax
        rhs(1,j) = ( - 25.d0 * fc(1,j) + 48.d0 * fc(2,j)
     &               - 36.d0 * fc(3,j) + 16.d0 * fc(4,j)
     &               -  3.d0 * fc(5,j) ) / ( 12.d0 * dx )
        
        rhs(2,j) = ( - 406.d0 * fc(1,j) - 300.d0 * fc(2,j)
     &               + 760.d0 * fc(3,j) -  80.d0 * fc(4,j)
     &               +  30.d0 * fc(5,j) -   4.d0 * fc(6,j) )
     &           /   ( 120.d0 * dx )
        
        do i = 3, ptsx-2
          rhs(i,j) = (           fc(i+2,j) - fc(i-2,j)
     &               + 28.d0 * ( fc(i+1,j) - fc(i-1,j) ) )
     &             / ( 12.d0 * dx )
        end do
        
        rhs(ptsx-1,j) = ( - 406.d0 * fc(ptsx,j)
     &                    - 300.d0 * fc(ptsx-1,j)
     &                    + 760.d0 * fc(ptsx-2,j)
     &                    -  80.d0 * fc(ptsx-3,j)
     &                    +  30.d0 * fc(ptsx-4,j)
     &                    -   4.d0 * fc(ptsx-5,j) )
     &                / ( - 120.d0 * dx )
        
        rhs(ptsx,j) = ( - 74.d0 * fc(ptsx,j)
     &                  + 16.d0 * fc(ptsx-1,j)
     &                  + 72.d0 * fc(ptsx-2,j)
     &                  - 16.d0 * fc(ptsx-3,j)
     &                  +  2.d0 * fc(ptsx-4,j) )
     &              / ( - 24.d0 * dx )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsxx(fc, rhs)
      
      ! RHS for the second derivative calculation in x direction
      implicit none
      include '../par.for'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do j = 1, jmax
        rhs(1,j) = (  9775.d0 * fc(1,j) - 20285.d0 * fc(2,j)
     &             + 11170.d0 * fc(3,j) -   550.d0 * fc(4,j)
     &             -   145.d0 * fc(5,j) +    35.d0 * fc(6,j) )
     &           / (    60.d0 * dxx )

        rhs(2,j) = ( 4834.d0 * fc(1,j) - 8424.d0 * fc(2,j)
     &             + 1890.d0 * fc(3,j) + 2320.d0 * fc(4,j)
     &             -  810.d0 * fc(5,j) +  216.d0 * fc(6,j)
     &             -   26.d0 * fc(7,j) ) / ( 360.d0 * dxx )

        do i = 3, ptsx-2
          rhs(i,j) = (   3.d0 * ( fc(i+2,j) + fc(i-2,j) )
     &               +  48.d0 * ( fc(i+1,j) + fc(i-1,j) )
     &               - 102.d0 *   fc(i,j) ) / ( 4.d0 * dxx )
        end do

       rhs(ptsx-1,j) = ( 4834.d0 * fc(ptsx,j)
     &                 - 8424.d0 * fc(ptsx-1,j)
     &                 + 1890.d0 * fc(ptsx-2,j)
     &                 + 2320.d0 * fc(ptsx-3,j)
     &                 -  810.d0 * fc(ptsx-4,j)
     &                 +  216.d0 * fc(ptsx-5,j)
     &                 -   26.d0 * fc(ptsx-6,j) )
     &               / (  360.d0 * dxx )
 
        rhs(ptsx,j) = (  9775.d0 * fc(ptsx,j)
     &                - 20285.d0 * fc(ptsx-1,j)
     &                + 11170.d0 * fc(ptsx-2,j)
     &                -   550.d0 * fc(ptsx-3,j)
     &                -   145.d0 * fc(ptsx-4,j)
     &                +    35.d0 * fc(ptsx-5,j) )
     &              / (    60.d0 * dxx )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsy(fc, rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do i = 1, ptsx
        rhs(i,1) = v_qdy(1) * ( fp_fd_coef_e(2) * fc(i,1)
     &                        + fp_fd_coef_e(3) * fc(i,2)
     &                        + fp_fd_coef_e(4) * fc(i,3)
     &                        + fp_fd_coef_e(5) * fc(i,4)
     &                        + fp_fd_coef_e(6) * fc(i,5) )
 
        rhs(i,2) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1)
     &                        + sp_fd_coef(5) * fc(i,2)
     &                        + sp_fd_coef(6) * fc(i,3)
     &                        + sp_fd_coef(7) * fc(i,4)
     &                        + sp_fd_coef(8) * fc(i,5)
     &                        + sp_fd_coef(9) * fc(i,6) )
 
       do j = 3, jmax-2
         rhs(i,j) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2)
     &                           + cp_fd_coef(5) * fc(i,j-1)
     &                           + cp_fd_coef(6) * fc(i,j)
     &                           + cp_fd_coef(7) * fc(i,j+1)
     &                           + cp_fd_coef(8) * fc(i,j+2) )
       end do
 
       rhs(i,jmax-1) = v_qdy(jmax-5) * ( pp_fd_coef(4)*fc(i,jmax)
     &                                 + pp_fd_coef(5)*fc(i,jmax-1)
     &                                 + pp_fd_coef(6)*fc(i,jmax-2)
     &                                 + pp_fd_coef(7)*fc(i,jmax-3)
     &                                 + pp_fd_coef(8)*fc(i,jmax-4)
     &                                 + pp_fd_coef(9)*fc(i,jmax-5) )
 
       rhs(i,jmax) = v_qdy(jmax-4) * ( lp_fd_coef(3) * fc(i,jmax)
     &                               + lp_fd_coef(4) * fc(i,jmax-1)
     &                               + lp_fd_coef(5) * fc(i,jmax-2)
     &                               + lp_fd_coef(6) * fc(i,jmax-3)
     &                               + lp_fd_coef(7) * fc(i,jmax-4) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine rhsyfv(fc, rhs)

      ! RHS for the first derivative calculation in y direction for uy
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do i = 1, ptsx

        rhs(i,1) = 0.d0
 
        rhs(i,2) = v_qdy(1) * ( sp_fd_coef(4) * fc(i,1)
     &                        + sp_fd_coef(5) * fc(i,2)
     &                        + sp_fd_coef(6) * fc(i,3)
     &                        + sp_fd_coef(7) * fc(i,4)
     &                        + sp_fd_coef(8) * fc(i,5)
     &                        + sp_fd_coef(9) * fc(i,6) )
 
       do j = 3, jmax-2
         rhs(i,j) = v_qdy(j-2) * ( cp_fd_coef(4) * fc(i,j-2)
     &                           + cp_fd_coef(5) * fc(i,j-1)
     &                           + cp_fd_coef(6) * fc(i,j)
     &                           + cp_fd_coef(7) * fc(i,j+1)
     &                           + cp_fd_coef(8) * fc(i,j+2) )
       end do
 
       rhs(i,jmax-1) = v_qdy(jmax-5) * ( pp_fd_coef(4)*fc(i,jmax)
     &                                 + pp_fd_coef(5)*fc(i,jmax-1)
     &                                 + pp_fd_coef(6)*fc(i,jmax-2)
     &                                 + pp_fd_coef(7)*fc(i,jmax-3)
     &                                 + pp_fd_coef(8)*fc(i,jmax-4)
     &                                 + pp_fd_coef(9)*fc(i,jmax-5) )
 
       rhs(i,jmax) = v_qdy(jmax-4) * ( lp_fd_coef(3) * fc(i,jmax)
     &                               + lp_fd_coef(4) * fc(i,jmax-1)
     &                               + lp_fd_coef(5) * fc(i,jmax-2)
     &                               + lp_fd_coef(6) * fc(i,jmax-3)
     &                               + lp_fd_coef(7) * fc(i,jmax-4) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyy(fc, rhs)

      ! RHS for the second derivative calculation in y direction
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j
      real*8 rhs(ptsx,jmax), fc(ptsx,jmax)

      do i = 1, ptsx

        rhs(i,1) = v_qdy2(1,1) * ( fp_sd_coef(3) * fc(i,1)
     &                           + fp_sd_coef(4) * fc(i,2)
     &                           + fp_sd_coef(5) * fc(i,3)
     &                           + fp_sd_coef(6) * fc(i,4)
     &                           + fp_sd_coef(7) * fc(i,5)
     &                           + fp_sd_coef(8) * fc(i,6) )

        rhs(i,2) = v_qdy2(1,1) * ( sp_sd_coef(4)  * fc(i,1)
     &                           + sp_sd_coef(5)  * fc(i,2)
     &                           + sp_sd_coef(6)  * fc(i,3)
     &                           + sp_sd_coef(7)  * fc(i,4)
     &                           + sp_sd_coef(8)  * fc(i,5)
     &                           + sp_sd_coef(9)  * fc(i,6)
     &                           + sp_sd_coef(10) * fc(i,7) )

       do j = 3, jmax-2
         rhs(i,j) = v_qdy2(j-2,1) * ( cp_sd_coef(4) * fc(i,j-2)
     &                              + cp_sd_coef(5) * fc(i,j-1)
     &                              + cp_sd_coef(6) * fc(i,j)
     &                              + cp_sd_coef(7) * fc(i,j+1)
     &                              + cp_sd_coef(8) * fc(i,j+2) )
       end do
 
      rhs(i,jmax-1) = v_qdy2(jmax-6,1)
     &              * ( pp_sd_coef(4)  * fc(i,jmax)
     &                + pp_sd_coef(5)  * fc(i,jmax-1)
     &                + pp_sd_coef(6)  * fc(i,jmax-2)
     &                + pp_sd_coef(7)  * fc(i,jmax-3)
     &                + pp_sd_coef(8)  * fc(i,jmax-4)
     &                + pp_sd_coef(9)  * fc(i,jmax-5)
     &                + pp_sd_coef(10) * fc(i,jmax-6) )
 
      rhs(i,jmax) = v_qdy2(jmax-5,1) * ( lp_sd_coef(3) * fc(i,jmax)
     &                                 + lp_sd_coef(4) * fc(i,jmax-1)
     &                                 + lp_sd_coef(5) * fc(i,jmax-2)
     &                                 + lp_sd_coef(6) * fc(i,jmax-3)
     &                                 + lp_sd_coef(7) * fc(i,jmax-4)
     &                                 + lp_sd_coef(8) * fc(i,jmax-5) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparxue(a, b, c, rhs)
 
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx), aux(2*inter), rhs(ptsx),
     &       u(ptsx), gam(ptsx), bet
 
      if (my_rank .eq. 0) then
        bet    = b(1)
        u(1)   = rhs(1) / bet
        gam(2) = c(1) / bet
        bet    = b(2) - a(2) * gam(2)
        u(2)   = ( rhs(2) - a(2) * u(1) ) / bet
       else
        call MPI_Recv(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank-1,
     &                57, MPI_COMM_WORLD, status, ierr)
        do i = 1, inter
          u(i)   = aux(i)
          gam(i) = aux(i+inter)
        end do
      end if
      i_ini = inter
      if (my_rank .eq. 0) i_ini = 3
      i_fim = ptsx - 2
      if (my_rank .eq. numproc) i_fim = ptsx
 
      bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
      do i = i_ini, i_fim
        gam(i) = c(i-1) / bet
        bet    = b(i) - a(i) * gam(i)
        u(i)   = ( rhs(i) - a(i) * u(i-1) ) / bet
      end do
 
      if (my_rank .lt. numproc) then
        do i = 1, inter
          aux(i)       = u(ptsx-inter+i-1)
          aux(i+inter) = gam(ptsx-inter+i-1)
        end do
        call MPI_Send(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank+1,
     &                57, MPI_COMM_WORLD, ierr)
      end if

      if (my_rank .lt. numproc) then
        call MPI_Recv(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank+1,
     &                67, MPI_COMM_WORLD, status, ierr)
        do i = 1, inter
          u(ptsx-inter+i) = aux(i)
          gam(ptsx-inter+i) = aux(i+inter)
        end do
      end if
 
      i_fim = ptsx - inter
      if (my_rank .eq. numproc) i_fim = ptsx - 1
      do i = i_fim, 1, -1
        u(i) = u(i) - gam(i+1) * u(i+1)
      end do

      if (my_rank .gt. 0) then
        do i = 1, inter
          aux(i)       = u(i+1)
          aux(i+inter) = gam(i+1)
        end do
        call MPI_Send(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank-1,
     &                67, MPI_COMM_WORLD, ierr)
      end if
      do i = 1, ptsx
        rhs(i) = u(i)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridparx(a, b, c, rhs)
 
      ! solves tridiagonal matrix for the derivatives in x direction
      ! domain decomposition in x direction -> parallel subroutine
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, i_ini, i_fim
      real*8 a(ptsx), b(ptsx), c(ptsx), aux(2*inter), rhs(ptsx,jmax),
     &       u(ptsx,jmax), gam(ptsx), bet
 
      do j = 1, jmax
        if (my_rank .eq. 0) then
          bet    = b(1)
          u(1,j) = rhs(1,j) / bet
          gam(2) = c(1) / bet
          bet    = b(2) - a(2) * gam(2)
          u(2,j) = ( rhs(2,j) - a(2) * u(1,j) ) / bet
         else
          call MPI_Recv(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank-1,
     &                  57, MPI_COMM_WORLD, status, ierr)
          do i = 1, inter
            u(i,j) = aux(i)
            gam(i) = aux(i+inter)
          end do
        end if
        i_ini = inter
        if (my_rank .eq. 0) i_ini = 3
        i_fim = ptsx - 2
        if (my_rank .eq. numproc) i_fim = ptsx
 
        bet = b(i_ini-1) - a(i_ini-1) * gam(i_ini-1)
        do i = i_ini, i_fim
          gam(i) = c(i-1) / bet
          bet    = b(i) - a(i) * gam(i)
          u(i,j) = ( rhs(i,j) - a(i) * u(i-1,j) ) / bet
        end do
 
        if (my_rank .lt. numproc) then
          do i = 1, inter
            aux(i)       = u(ptsx-inter+i-1,j)
            aux(i+inter) = gam(ptsx-inter+i-1)
          end do
          call MPI_Send(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank+1,
     &                  57, MPI_COMM_WORLD, ierr)
        end if
      end do

      do j = 1, jmax
        if (my_rank .lt. numproc) then
          call MPI_Recv(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank+1,
     &                  67, MPI_COMM_WORLD, status, ierr)
          do i = 1, inter
            u(ptsx-inter+i,j) = aux(i)
            gam(ptsx-inter+i) = aux(i+inter)
          end do
        end if
 
        i_fim = ptsx - inter
        if (my_rank .eq. numproc) i_fim = ptsx - 1
        do i = i_fim, 1, -1
          u(i,j) = u(i,j) - gam(i+1) * u(i+1,j)
        end do

        if (my_rank .gt. 0) then
          do i = 1, inter
            aux(i)       = u(i+1,j)
            aux(i+inter) = gam(i+1)
          end do
          call MPI_Send(aux, 2*inter, MPI_DOUBLE_PRECISION, my_rank-1,
     &                  67, MPI_COMM_WORLD, ierr)
        end if
        do i = 1, ptsx
          rhs(i,j) = u(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridy(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include '../par.for'
      integer i, j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet,
     &       rhs(ptsx,jmax), u(jmax)

      do i = 1, ptsx
        bet  = b(1)
        u(1) = rhs(i,1) / bet
        do j = 2, jmax
          gam(j) = c(j-1) / bet
          bet    = b(j) - a(j) * gam(j)
          u(j)   = ( rhs(i,j) - a(j) * u(j-1) ) / bet
        end do
        do j = jmax - 1, 1, -1
          u(j) = u(j) - gam(j+1) * u(j+1)
        end do
        do j = 1, jmax
          rhs(i,j) = u(j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridseqx(a, b, c, rhs)

      ! solves tridiagonal matrix for the derivatives in x direction
      implicit none
      include '../par.for'
      integer i, j
      real*8 a(ptsx), b(ptsx), c(ptsx), gam(ptsx), bet,
     &       rhs(ptsx,jmax), u(ptsx)

      do j = 1, jmax
        bet  = b(1)
        u(1) = rhs(1,j) / bet
        do i = 2, ptsx
          gam(i) = c(i-1) / bet
          bet    = b(i) - a(i) * gam(i)
          u(i)   = ( rhs(i,j) - a(i) * u(i-1) ) / bet
        end do
        do i = ptsx - 1, 1, -1
          u(i) = u(i) - gam(i+1) * u(i+1)
        end do
        do i = 1, ptsx
          rhs(i,j) = u(i)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coef(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_fd_coef_e(1)
      c(1)      = 0.d0

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax-2
        a(j)    = cp_fd_coef(1)
        b(j)    = cp_fd_coef(2)
        c(j)    = cp_fd_coef(3)
      end do

      a(jmax-1) = pp_fd_coef(3)
      b(jmax-1) = pp_fd_coef(2)
      c(jmax-1) = pp_fd_coef(1)

      a(jmax)   = lp_fd_coef(2)
      b(jmax)   = lp_fd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeffv(a, b, c)

      ! mount the LHS of the matrix for the first derivative of uy
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax-2
        a(j)    = cp_fd_coef(1)
        b(j)    = cp_fd_coef(2)
        c(j)    = cp_fd_coef(3)
      end do

      a(jmax-1) = pp_fd_coef(3)
      b(jmax-1) = pp_fd_coef(2)
      c(jmax-1) = pp_fd_coef(1)

      a(jmax)   = lp_fd_coef(2)
      b(jmax)   = lp_fd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeff(a, b, c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_sd_coef(1)
      c(1)      = fp_sd_coef(2)

      a(2)      = sp_sd_coef(1)
      b(2)      = sp_sd_coef(2)
      c(2)      = sp_sd_coef(3)

      do j = 3, jmax-2
        a(j)    = cp_sd_coef(1)
        b(j)    = cp_sd_coef(2)
        c(j)    = cp_sd_coef(3)
      end do

      a(jmax-1) = pp_sd_coef(3)
      b(jmax-1) = pp_sd_coef(2)
      c(jmax-1) = pp_sd_coef(1)

      a(jmax)   = lp_sd_coef(2)
      b(jmax)   = lp_sd_coef(1)
      c(jmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coefx(a, b, c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include '../par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, ptsx-2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(ptsx-1) = 2.d0
      b(ptsx-1) = 6.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 4.d0
      b(ptsx)   = 1.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeffx(a, b, c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include '../par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 13.d0
      c(1)      = 137.d0

      a(2)      = 1.d0
      b(2)      = 12.d0
      c(2)      = 3.d0

      do i = 3, ptsx-2
        a(i)    = 2.d0
        b(i)    = 11.d0
        c(i)    = 2.d0
      end do

      a(ptsx-1) = 3.d0
      b(ptsx-1) = 12.d0
      c(ptsx-1) = 1.d0

      a(ptsx)   = 137.d0
      b(ptsx)   = 13.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_derivs(var)

      ! exchange values of the boundaries
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, j
      real*8 aux(meshdx*jmax), var(ptsx,jmax)

      ! variable used to calculate the number of columns needed
      ! to go forward or backward (interm)

      if (my_rank .lt. numproc) then
        ! Sending the near right boundary columns to node + 1
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax) = var(i+ptsx-24-1,j)
          end do
        end do
        call MPI_Send(aux, meshdx*jmax, MPI_DOUBLE_PRECISION,
     &                my_rank+1, 10, MPI_COMM_WORLD, ierr)
        ! Receiving the new right boundary columns from node + 1
        call MPI_Recv(aux, meshdx*jmax, MPI_DOUBLE_PRECISION,
     &                my_rank+1, 20, MPI_COMM_WORLD, status, ierr)
        do i = 1, meshdx
          do j = 1, jmax
            var(i+ptsx-meshdx,j) = aux(j+(i-1)*jmax)
          end do
        end do
 
       end if

      if (my_rank .gt. 0) then
        ! Receiving the new left boundary columns from node - 1
        call MPI_Recv(aux, meshdx*jmax, MPI_DOUBLE_PRECISION,
     &                my_rank-1, 10, MPI_COMM_WORLD, status, ierr)
        do i = 1, meshdx
          do j = 1, jmax
            var(i,j) = aux(j+(i-1)*jmax)
          end do
        end do
        ! Sending the new right boundary columns to node - 1
        do i = 1, meshdx
          do j = 1, jmax
            aux(j+(i-1)*jmax) = var(i+24-meshdx+1,j)
          end do
        end do
        call MPI_Send(aux, meshdx*jmax, MPI_DOUBLE_PRECISION,
     &                my_rank-1, 20, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         end of derivative calculations                       c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
