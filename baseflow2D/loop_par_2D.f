cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                            loop calculations                                 c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine outuy

      ! outflow uy calculation with 6th order compact appx.
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      integer j, indx(jmax)
      real*8 a(jmax,5), al(jmax,5), dwzdx(ptsx,jmax), rhs(jmax), dya
      common/dwdx/ dwzdx

      call derparx(dwzdx, wz)

      if (my_rank .eq. numproc) then

        call lhsoutuy(a)
        call bandy5(a, al, indx)
        rhs(1) = 0.d0
        rhs(2) = - dy * dy * ( sp_poi_coef(2,1) * dwzdx(ptsx,2)
     &                       + sp_poi_coef(3,1) * dwzdx(ptsx,3) )
        do j = 3, jmax - 2
          dya = dy * stf**(j-3)
          rhs(j) = - dya*dya * ( cp_poi_coef(1,1) * dwzdx(ptsx,j-1)
     &                         + cp_poi_coef(2,1) * dwzdx(ptsx,j)
     &                         + cp_poi_coef(3,1) * dwzdx(ptsx,j+1) )
        end do
        dya          =  dy * (stf**(jmax-3))
        rhs(jmax-1) = - dya*dya * pp_poi_coef(1,1)*dwzdx(ptsx,jmax-1)
        rhs(jmax)   = - dya*dya * lp_poi_coef(1,1)*dwzdx(ptsx,jmax)
     &                + dya * lp_poi_coef(2,1) * duexmdx(ptsx)
      
c       call banbky5(a, al, indx, rhs)
        do j = 2, jmax
c         uy(ptsx,j) = rhs(j)
          uy(ptsx,j) = ( 770.d0 * uy(ptsx-1,j) - 1070.d0 * uy(ptsx-2,j)
     &                 + 780.d0 * uy(ptsx-3,j) -  305.d0 * uy(ptsx-4,j)
     &                 +  50.d0 * uy(ptsx-5,j) ) / 225.d0
        end do

      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lhsoutuy(a)

      ! LHS of the uy calculation
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer j
      real*8 a(jmax,5)

      a(1,1)      = 0.d0
      a(1,2)      = 0.d0
      a(1,3)      = 1.d0
      a(1,4)      = 0.d0
      a(1,5)      = 0.d0

      a(2,1)      = 0.d0
      a(2,2)      = sp_poi_coef(4,1)
      a(2,3)      = sp_poi_coef(5,1)
      a(2,4)      = sp_poi_coef(6,1)
      a(2,5)      = sp_poi_coef(7,1)

      do j = 3, jmax - 2
        a(j,1)    = cp_poi_coef(4,1)
        a(j,2)    = cp_poi_coef(5,1)
        a(j,3)    = cp_poi_coef(6,1)
        a(j,4)    = cp_poi_coef(7,1)
        a(j,5)    = cp_poi_coef(8,1)
      end do

      a(jmax-1,1) = 0.d0
      a(jmax-1,2) = pp_poi_coef(2,1)
      a(jmax-1,3) = pp_poi_coef(3,1)
      a(jmax-1,4) = pp_poi_coef(4,1)
      a(jmax-1,5) = 0.d0

      a(jmax,1)   = lp_poi_coef(5,1)
      a(jmax,2)   = lp_poi_coef(4,1)
      a(jmax,3)   = lp_poi_coef(3,1)
      a(jmax,4)   = 0.d0
      a(jmax,5)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine wz_wall

      ! calculate the vorticity in z direction at the wall using the
      ! v-poisson equation
      ! the variable lapv comes from subroutine wx_wall
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i
      real*8 rhs(ptsx), aux(2)

      do i = 1, ptsx
        rhs(i) = - ( w_at_w_coef(3) * uy(i,1)
     &             + w_at_w_coef(4) * uy(i,2)
     &             + w_at_w_coef(5) * uy(i,3)
     &             + w_at_w_coef(6) * uy(i,4)
     &             + w_at_w_coef(7) * uy(i,5)
     &             + w_at_w_coef(8) * uy(i,6)
     &             + w_at_w_coef(9) * uy(i,7) )
     &             / ( dyy * w_at_w_coef(1) )
      end do

      if (my_rank .gt. 0) then
        call MPI_Recv(aux, 2, mpi_double_precision, my_rank - 1,
     &                571, MPI_COMM_WORLD, status, ierr)
        wz(1,1) = aux(1)
        wz(2,1) = aux(2)
       else
        wz(2,1) = ( dx * ( 251.d0 * rhs(1)   + 646.d0 * rhs(2)
     &                   - 264.d0 * rhs(3)   + 106.d0 * rhs(4)
     &                   -  19.d0 * rhs(5) ) + 720.d0 * wz(1,1) )
     &            / 720.d0
      end if

      do i = 3, ptsx - 2
        wz(i,1) = ( dx * (  281.d0 * rhs(i-2) + 2056.d0 * rhs(i-1)
     &                   + 1176.d0 * rhs(i)   -  104.d0 * rhs(i+1)
     &                   +   11.d0 * rhs(i+2) ) / 90.d0
     &            + 11.d0 * wz(i-2,1) + 16.d0 * wz(i-1,1) )
     &            / 27.d0
      end do

      if (my_rank .lt. numproc) then
        aux(1) = wz(ptsx-inter,1)
        aux(2) = wz(ptsx-inter+1,1)
        call MPI_Send(aux, 2, mpi_double_precision, my_rank + 1,
     &                571, MPI_COMM_WORLD, ierr)
       else
        i = ptsx - 1
        wz(i,1) = ( dx * ( 10.d0 * rhs(i-2) + 57.d0 * rhs(i-1)
     &                   + 24.d0 * rhs(i)   -         rhs(i+1) )
     &            + 33.d0 * wz(i-2,1) + 24.d0 * wz(i-1,1) )
     &            / 57.d0
        i = ptsx
        wz(i,1) = ( dx * ( 251.d0 * rhs(i)   + 646.d0 * rhs(i-1)
     &                   - 264.d0 * rhs(i-2) + 106.d0 * rhs(i-3)
     &                   -  19.d0 * rhs(i-4) )
     &            + 720.d0 * wz(i-1,1) ) / 720.d0
       end if

!!!!!!!!! VER SE PRECISA DISTO MESMO
      call boundary_exchange_wwall(wz)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine poi_ux

      ! calculate the new ux velocity using the continuity equation
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include '../comm.fs'
      include 'mpif.h'
      integer status(MPI_status_size)
      integer i, j, i_ini
      real*8 var(jmax), xad, ue, aux(2), dvdx(ptsx,jmax),
     &       duydy(ptsx,jmax), dvdxa, dya, m, ueptsx(ptsx)

!     if (my_rank .eq. 0) then
!       do j = 1, jmax
!         dvdxa = ( - 1764.d0 * uy(1,j) + 4320.d0 * uy(2,j)
!    &              - 5400.d0 * uy(3,j) + 4800.d0 * uy(4,j)
!    &              - 2700.d0 * uy(5,j) + 864.d0 *  uy(6,j)
!    &              -  120.d0 * uy(7,j) ) / (720.d0 * dx)
!         var(j) = dvdxa / fac_y  + wz(1,j)
!       end do
!    
!       dya = dy
!       ux(1,2) = ( dya * ( sp_integ_coef(3) * var(1)
!    &                    + sp_integ_coef(4) * var(2)
!    &                    + sp_integ_coef(5) * var(3)
!    &                    + sp_integ_coef(6) * var(4)
!    &                    + sp_integ_coef(7) * var(5)  )
!    &                    + sp_integ_coef(1) * ux(1,1) )
!    &            / ( - 1.d0 * sp_integ_coef(2) )
!       do j = 3, jmax - 2
!          dya = dy * stf**(j-3)
!          ux(1,j) = ( dya * ( cp_integ_coef(4) * var(j-2)
!    &                       + cp_integ_coef(5) * var(j-1)
!    &                       + cp_integ_coef(6) * var(j)
!    &                       + cp_integ_coef(7) * var(j+1)
!    &                       + cp_integ_coef(8) * var(j+2)  )
!    &                       + cp_integ_coef(1) * ux(1,j-2)
!    &                       + cp_integ_coef(2) * ux(1,j-1) )
!    &               / ( - 1.d0 * cp_integ_coef(3) )
!       end do
!
!       j = jmax - 1
!       dya = dy * stf**(j-3)
!       ux(1,j) = ( dya * ( pp_integ_coef(4) * var(j-2)
!    &                    + pp_integ_coef(5) * var(j-1)
!    &                    + pp_integ_coef(6) * var(j)
!    &                    + pp_integ_coef(7) * var(j+1)  )
!    &                    + pp_integ_coef(1) * ux(1,j-2)
!    &                    + pp_integ_coef(2) * ux(1,j-1) )
!    &            / ( - 1.d0 * pp_integ_coef(3) )

!       dya = dy * stf**(jmax-5)
!       ux(1,jmax) = ( dya * ( lp_integ_coef(3) * var(jmax)
!    &                       + lp_integ_coef(4) * var(jmax-1)
!    &                       + lp_integ_coef(5) * var(jmax-2)
!    &                       + lp_integ_coef(6) * var(jmax-3)
!    &                       + lp_integ_coef(7) * var(jmax-4)  )
!    &                       + lp_integ_coef(2) * ux(1,jmax-1) )
!    &               / ( - 1.d0 * lp_integ_coef(1) )
!
!     end if

      call deryfv(duydy, uy)
      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      do j = 2, jmax - 1
        if (my_rank .gt. 0) then
          call MPI_Recv(aux, 2, mpi_double_precision, my_rank - 1,
     &                  153, MPI_COMM_WORLD, status, ierr)
          ux(1,j) = aux(1)
          ux(2,j) = aux(2)
         else
          ux(2,j) = ( - dx * ( 251.d0 * duydy(1,j)
     &                       + 646.d0 * duydy(2,j)
     &                       - 264.d0 * duydy(3,j)
     &                       + 106.d0 * duydy(4,j)
     &                       -  19.d0 * duydy(5,j) )
     &                       + 720.d0 * ux(1,j) ) / 720.d0
        end if

        do i = 3, ptsx - 2
          ux(i,j) = ( - dx * (  281.d0 * duydy(i-2,j)
     &                       + 2056.d0 * duydy(i-1,j)
     &                       + 1176.d0 * duydy(i,j)
     &                       -  104.d0 * duydy(i+1,j)
     &                       +   11.d0 * duydy(i+2,j) ) / 90.d0
     &                + 11.d0 * ux(i-2,j) + 16.d0 * ux(i-1,j) )
     &                / 27.d0
        end do

        if (my_rank .lt. numproc) then
          aux(1) = ux(ptsx-inter,j)
          aux(2) = ux(ptsx-inter+1,j)
          call MPI_Send(aux, 2, mpi_double_precision, my_rank + 1,
     &                  153, MPI_COMM_WORLD, ierr)
         else
          i = ptsx - 1
          ux(i,j) = ( - dx * ( 10.d0 * duydy(i-2,j)
     &                       + 57.d0 * duydy(i-1,j)
     &                       + 24.d0 * duydy(i,j)
     &                       -         duydy(i+1,j) )
     &                + 33.d0 * ux(i-2,j) + 24.d0 * ux(i-1,j) )
     &                / 57.d0
          i = ptsx
          ux(i,j) = ( - dx * ( 251.d0 * duydy(i,j)
     &                       + 646.d0 * duydy(i-1,j)
     &                       - 264.d0 * duydy(i-2,j)
     &                       + 106.d0 * duydy(i-3,j)
     &                       -  19.d0 * duydy(i-4,j) )
     &                       + 720.d0 * ux(i-1,j) ) / 720.d0
        end if
      end do

!!!!!!!!! VER SE PRECISA DISTO MESMO
      call boundary_exchange_1stmode(ux)

!     if (my_rank .eq. 0) then
!       ue = ux(1,jmax)
!     end if
!     call MPI_BCAST(ue, 1, mpi_double_precision, 0, mpi_comm_world,
!    &               ierr)
!     do i = 1, ptsx
!       xad        = dble(i+shift-1) * dx + x0
!       m          = beta_fs(i+shift) / (2.d0 - beta_fs(i+shift))
!       ux(i,jmax) = ue * xad**m
!       ueptsx(i)  = ux(i,jmax)
!     end do
!     call derparxue(duexmdx, ueptsx)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_wwall(var)

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer i
      integer status(MPI_Status_size)
      real*8 aux(3), var(ptsx,jmax)

      if (my_rank .lt. numproc) then
        ! receives the lasts columns
        call MPI_Recv(aux, 3, mpi_double_precision, my_rank + 1,
     &                140, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          var(i + ptsx - 3, 1) = aux(i)
        end do
      end if
      if (my_rank .gt. 0) then
        ! sends the seconds columns
        do i = 1, 3
          aux(i) = var(i + inter - 2, 1)
        end do
        call MPI_Send(aux, 3, mpi_double_precision, my_rank - 1,
     &                140, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_1stmode(var)

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer i, j
      integer status(MPI_Status_size)
      real*8 aux(3*jmax), var(ptsx,jmax)

      if (my_rank .lt. numproc) then
        ! receives the lasts columns
        call MPI_Recv(aux, 3*jmax, mpi_double_precision, my_rank + 1,
     &                240, MPI_COMM_WORLD, status, ierr)
        do i = 1, 3
          do j = 1, jmax
            var(i + ptsx - 3,j) = aux(j+(i-1)*jmax)
          end do
        end do
      end if
      if (my_rank .gt. 0) then
        ! sends the seconds columns
        do i = 1, 3
          do j = 1, jmax
            aux(j+(i-1)*jmax) = var(i + inter - 2,j)
          end do
        end do
        call MPI_Send(aux, 3*jmax, mpi_double_precision, my_rank - 1,
     &                240, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine filter_trid

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      integer i, j, i_ini
      real*8 afil(ptsx), bfil(ptsx), cfil(ptsx), rhsz(ptsx,jmax)
      common/filt/ afil, bfil, cfil

      i_ini = 1
      if (my_rank.eq.0) i_ini = 2

      call rhsf_trid(rhsz)
      call tridseqx(afil, bfil, cfil, rhsz)
      call boundary_exchange_derivs(rhsz)
      do i = i_ini, ptsx
        do j = 2, jmax
          wz(i,j) = rhsz(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lhs_tridf(a, b, c)

      ! mount the LHS of the matrix for the second derivative
      implicit none
      include '../par.for'
      integer i
      real*8 a(ptsx), b(ptsx), c(ptsx)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 0.d0

      a(2)      = 0.d0
      b(2)      = 1.d0
      c(2)      = 0.d0

      a(3)      = 0.d0
      b(3)      = 1.d0
      c(3)      = 0.d0

      do i = 4, ptsx - 3
        a(i)    = alphaf
        b(i)    = 1.d0
        c(i)    = alphaf
      end do

      a(ptsx-2) = 0.d0
      b(ptsx-2) = 1.d0
      c(ptsx-2) = 0.d0

      a(ptsx-1) = 0.d0
      b(ptsx-1) = 1.d0
      c(ptsx-1) = 0.d0

      a(ptsx)   = 0.d0
      b(ptsx)   = 1.d0
      c(ptsx)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsf_trid(rhsz)

      ! calculate the RHS for the wx filter
      implicit none
      include '../par.for'
      include 'comm.var'
      integer i, j
      real*8 rhsz(ptsx,jmax)

      do j = 1, jmax
        rhsz(1,j) = (   15.d0 * wz(1,j) +  4.d0 * wz(2,j)
     &                -  6.d0 * wz(3,j) +  4.d0 * wz(4,j)
     &                -         wz(5,j) ) / 16.d0
        rhsz(2,j) = (           wz(1,j) + 12.d0 * wz(2,j)
     &                +  6.d0 * wz(3,j) -  4.d0 * wz(4,j)
     &                +         wz(5,j) ) / 16.d0
        rhsz(3,j) = ( -         wz(1,j) +  4.d0 * wz(2,j)
     &                + 10.d0 * wz(3,j) +  4.d0 * wz(4,j)
     &                -         wz(5,j) ) / 16.d0

        do i = 4, ptsx - 3
          rhsz(i,j) = af *   wz(i,j)
     &              + bf * ( wz(i+1,j) + wz(i-1,j))
     &              + cf * ( wz(i+2,j) + wz(i-2,j))
     &              + df * ( wz(i+3,j) + wz(i-3,j))
        end do

        rhsz(ptsx-2,j) = ( -         wz(ptsx,j)
     &                     +  4.d0 * wz(ptsx-1,j)
     &                     + 10.d0 * wz(ptsx-2,j)
     &                     +  4.d0 * wz(ptsx-3,j)
     &                     -         wz(ptsx-4,j) ) / 16.d0
        rhsz(ptsx-1,j) = (           wz(ptsx,j)
     &                     + 12.d0 * wz(ptsx-1,j)
     &                     +  6.d0 * wz(ptsx-2,j)
     &                     -  4.d0 * wz(ptsx-3,j)
     &                     +         wz(ptsx-4,j) ) / 16.d0
        rhsz(ptsx,j)   = (   15.d0 * wz(ptsx,j)
     &                     +  4.d0 * wz(ptsx-1,j)
     &                     -  6.d0 * wz(ptsx-2,j)
     &                     +  4.d0 * wz(ptsx-3,j)
     &                     -         wz(ptsx-4,j) ) / 16.d0

      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         end of loop calculations                             c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
