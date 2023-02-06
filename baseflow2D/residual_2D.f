
      subroutine ns_residual(t)

      ! NS-residual residual calculation

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include 'mpif.h'

      integer :: i, j, t, it
      integer :: status(MPI_status_size)
      double precision :: local_max(3), global_max(3)
      double precision ::      res(ptsx,jmax),      S_z(ptsx,jmax)
      double precision ::    dwzdx(ptsx,jmax),    dwzdy(ptsx,jmax)
      double precision ::  d2wzdx2(ptsx,jmax),  d2wzdy2(ptsx,jmax)
      double precision ::   duxdxe(ptsx,jmax),   duxdye(ptsx,jmax)
      double precision ::   duydxe(ptsx,jmax),   duydye(ptsx,jmax)
      double precision ::   dwzdxe(ptsx,jmax),   dwzdye(ptsx,jmax)
      double precision :: d2uxdx2e(ptsx,jmax), d2uxdy2e(ptsx,jmax)
      double precision :: d2uydx2e(ptsx,jmax), d2uydy2e(ptsx,jmax)
      double precision :: d2wzdx2e(ptsx,jmax), d2wzdy2e(ptsx,jmax)
      common/mmsterms/ S_z
      common/residual/ res
      common/dexact/ duxdxe, duxdye, duydxe, duydye, dwzdxe, dwzdye,
     &               d2uxdx2e, d2uxdy2e, d2uydx2e, d2uydy2e,
     &               d2wzdx2e, d2wzdy2e

      ! x-derivative calculation
      call derparx(dwzdx, wz)
      call derparxx(d2wzdx2, wz)

      ! y-derivative calculation
      call dery(dwzdy, wz)
      call deryy(d2wzdy2, wz)

      ! residual calculation
      local_max = 0.d0
      do j = 1, jmax
        do i = 1, ptsx
          res(i,j) = ux(i,j) * dwzdx(i,j) + uy(i,j) * dwzdy(i,j)
     &             - ( d2wzdx2(i,j) + fac_y * d2wzdy2(i,j) ) / Re
     &             - S_z(i,j)
           if (local_max(3) .lt. abs(res(i,j))) then
             local_max(1) = i
             local_max(2) = j
             local_max(3) = abs(res(i,j))
           end if
        end do
      end do

      if (my_rank .gt. 0) then
        ! Send the maximum residual to the node 1
        call MPI_Send(local_max, 3, MPI_DOUBLE_PRECISION, 0, 51,
     &                MPI_COMM_WORLD, ierr)
       else
        global_max = local_max
        do it = 1, numproc
          ! Receive residual from all processing elements and find the global maximum
          call MPI_Recv(local_max, 3, MPI_DOUBLE_PRECISION, it, 51,
     &                  MPI_COMM_WORLD, status, ierr)
          if (local_max(3) .gt. global_max(3)) then
            global_max = local_max
          end if
        end do
        write(*,*) 'residual', t, int(global_max(1)),
     &             int(global_max(2)), global_max(3)
        if (global_max(3) < 1d-9) then
          write(*,*) 'Stationary state.'
          call escreve
          call escreve_residual(t)
          return
        end if
      end if

      return
      end subroutine ns_residual

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve_residual(t)

      ! write residual

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include 'mpif.h'

      character*40 file_name
      integer :: i, j, t
      double precision :: x, y, res(ptsx,jmax)
      common/residual/ res

      write(file_name,'(a,i0.5,a,i0.2,a)') 'res_',t,'_',my_rank,'.dat'
      open (3, file = file_name, status = 'unknown')
      write(3,*) 'VARIABLES="x","y","residual"'
      write(3,*) 'ZONE I=',ptsx,', J=',jmax,', F=POINT'
      do j = 1, jmax
        if (stf .eq. 1.d0) then
          ! without stretching
          y = dble(j-1) * dy
         else
          ! with stretching
          y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
        end if
        y = y / sqrt(fac_y)
        do i = 1, ptsx
          x = x0 + dble(i-1+shift) * dx
          write(3,*) x, y, abs(res(i,j))
        end do
      end do
      close(3)

      return
      end subroutine escreve_residual
