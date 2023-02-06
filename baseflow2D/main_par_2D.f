cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                              main subroutines                                c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program baseflow_20140819

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      logical run
      integer p
      real*8 tempo_initial, tempo_final

      ! start up MPI
      call MPI_Init(ierr)

      ! find out process rank
      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

      ! find out number of processes
      call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)

      ! p nodes in the environment. So the number of the last node is numproc
      numproc = p - 1

      ! calculate the i_shift from one computer to another with the domain
      ! decomposition
      shift = my_rank * (ptsx - inter - 1)

      ! to calculate the total time simulation
      if (my_rank .eq. 0) then
        ! get the initial time
        tempo_initial = mpi_wtime()
      end if

      ! verification of adopted points in x and y directions
      run = .true.
      call verifica(run)
      if (run) call solver

      if (my_rank .eq. 0) then
        ! get the final time
        tempo_final = mpi_wtime()
        ! calculate the total time of simulation
        tempo_final = tempo_final - tempo_initial
        write(*,4) tempo_final
      end if

      ! finalize MPI
      call MPI_Finalize(ierr)

   4  format(1x, 1e21.9e3)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine verifica(run)

      implicit none
      include '../par.for'
      include '../comm.par'
      logical run
      integer nodesx, chute1, chute2
      real*8 n_npr

      ! calculate the number of points in the x direction for each node
      nodesx =     ( imax + (inter + 1) * numproc ) /     (numproc + 1)
      n_npr  = dble( imax + (inter + 1) * numproc ) / dble(numproc + 1)

      ! stop the program if the number of total points is not exactly
      ! divided by nodes
      if (ptsx .ne. nodesx) then
        write(*,*) 'ALTERAR O VALOR DE PTSX NO ARQUIVO par.for PARA:'
        write(*,*) nodesx
        run = .false.
      end if

      ! stop the program if the number of points in the x-direction can
      ! not be used in multigrid program 
      if ( (nodesx .ne. n_npr) .or.
     &     (mod(nodesx-1, 2**(msh-1)) .ne. 0) ) then
        chute1 = (   (nodesx - 1) / 2**(msh - 1) )       * 2**(msh - 1)
        chute2 = ( ( (nodesx - 1) / 2**(msh - 1) ) + 1 ) * 2**(msh - 1)
        write(*,*)' Number of points in the x direction can
     &              not be used in multigrid solver'
        write(*,*)' The number of points in the x direction
     &              should be:'
        write(*,*)  chute1 * (numproc + 1) - numproc * inter + 1,'   or'
        write(*,*)  chute2 * (numproc + 1) - numproc * inter + 1
        run = .false.
      end if

      ! stop the program if the number of points in the x-direction can
      ! not be used in multigrid program 
      if ( (numproc .gt. 1) .and. (nodesx .lt. 2*inter) ) then
        write(*,*) ' Number of points in the x direction can not be
     &               used in multigrid solver'
        write(*,*) ' The number of points in the x direction should
     &               be:'
        write(*,*) ((nodesx * 2) - 1) * (numproc + 1) - numproc
     &             * inter + 2
        run = .false.
      end if

      ! stop the program if the number of points in the y-direction can
      ! not be used in multigrid program
      if ( (dble(jmax-1) / dble(2**(msh-1))) .ne.
     &     ((jmax-1) / 2**(msh-1)) ) then
        write(*,*)' Number of points in the y direction can
     &              not be used in multigrid solver'
        write(*,*)' The number of points in the y direction
     &              should be ', 2**(msh-1) *
     &              ((jmax - 1) / 2**(msh-1)) + 1,' or ',
     &              2**(msh-1) * (((jmax - 1) / 2**(msh-1)) + 1) + 1
        run = .false.
      end if

      ! stop the program if the number of points in the y-direction can
      ! not be used in multigrid program
      if (jmax .lt. 2**msh) then
        write(*,*) ' Number of points in the y direction can not be
     & used in multigrid solver'
        write(*,*) ' The number of points in the y direction should
     & be:'
        write(*,*) jmax * 2 - 1
        run = .false.
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine solver

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include 'mpif.h'
      integer i, j, t, i_ini, seed
      real*8 dv1(ptsx,jmax), wz1(ptsx,jmax), dvt1(ptsx,jmax),
     &       th1(ptsx,jmax), ux_exact(ptsx,jmax), rd,
     &       uy_exact(ptsx,jmax), wz_exact(ptsx,jmax),
     &       maxlocal_wz, maxlocal_th
      common/exact/ ux_exact, uy_exact, wz_exact

      ! initial point of the domain for the solution computation
      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      select case (my_form)

        case(1) ! -> Gortler vortices simulations without heat transfer

          ! initialize variables
          call init(dv1, wz1)
          ! calculate the initial residual
c          call ns_residual(0)

          ! time integration
          do t = 1, tt_base
            ! fisrt step of second-order Runge-Kutta
            wz(:,jmax) = 0.0d0
            wz1 = wz
            call drv(dv1)
c            dv1_old = dv1
            do j = 2, jmax-1
              do i = i_ini, ptsx
                wz(i,j)  = wz1(i,j) + 0.5*dt_base * dv1(i,j)
              end do
              if (my_rank.eq.numproc) then
          wz(ptsx,j) = ( 770.d0 * wz(ptsx-1,j) - 1070.d0 * wz(ptsx-2,j)
     &                 + 780.d0 * wz(ptsx-3,j) -  305.d0 * wz(ptsx-4,j)
     &                 +  50.d0 * wz(ptsx-5,j) ) / 225.d0
              end if
            end do
            call loop(1d-6)

            ! second step of second-order Runge-Kutta
            call drv(dv1)
            do j = 2, jmax-1
              do i = i_ini, ptsx
               wz(i,j) = wz1(i,j) + dt_base*dv1(i,j)
              end do
              if (my_rank.eq.numproc) then
          wz(ptsx,j) = ( 770.d0 * wz(ptsx-1,j) - 1070.d0 * wz(ptsx-2,j)
     &                 + 780.d0 * wz(ptsx-3,j) -  305.d0 * wz(ptsx-4,j)
     &                 +  50.d0 * wz(ptsx-5,j) ) / 225.d0
              end if
            end do
            call loop(1d-7)

            write(*,*) my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax)
            if (mod(t, 5000) .eq. 0) call escreve(t)

            if (mod(t,100) .eq. 0) then
              maxlocal_wz = maxval(dabs(wz1(:,:) - wz(:,:)))
              write(*,*) my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax),
     &                   maxlocal_wz
              call isstat(maxlocal_wz)
              if ( (maxlocal_wz .lt. 1d-8) .and. (t .ge. stpp) ) then
                write(*,*) 'Stationary state.'
                call escreve2(t)
                return
              end if
            end if

c            write(*,*) my_rank, t            
          end do

          call escreve(t)
 
        case(2) ! -> Gortler vortices simulations with heat transfer 
 
         ! initialize variables
          call init_theta(dv1, wz1, dvt1, th1)

          ! time integration
          do t = 1, tt_base
            wz(:,jmax) = 0.0d0
            th(:,jmax) = 1.0d0
            wz1 = wz
            th1 = th
            ! first step of second-order Runge-Kutta
            call drv_theta(dv1, dvt1)
            do j = 2, jmax-1
              do i = i_ini, ptsx
                wz(i,j)  = wz1(i,j) + 0.5*dt_base * dv1(i,j)
                th(i,j)  = th1(i,j) + 0.5*dt_base * dvt1(i,j)
              end do
              if (my_rank.eq.numproc) then
          wz(ptsx,j) = ( 770.d0 * wz(ptsx-1,j) - 1070.d0 * wz(ptsx-2,j)
     &                 + 780.d0 * wz(ptsx-3,j) -  305.d0 * wz(ptsx-4,j)
     &                 +  50.d0 * wz(ptsx-5,j) ) / 225.d0
          th(ptsx,j) = ( 770.d0 * th(ptsx-1,j) - 1070.d0 * th(ptsx-2,j)
     &                 + 780.d0 * th(ptsx-3,j) -  305.d0 * th(ptsx-4,j)
     &                 +  50.d0 * th(ptsx-5,j) ) / 225.d0
              end if
            end do
            call loop(1d-6)

            ! second step of second-order Runge-Kutta
            call drv_theta(dv1, dvt1)
            do j = 2, jmax
              do i = i_ini, ptsx
                wz(i,j)  = wz1(i,j) + dt_base * dv1(i,j)
                th(i,j)  = th1(i,j) + dt_base * dvt1(i,j)
              end do
              if (my_rank.eq.numproc) then
          wz(ptsx,j) = ( 770.d0 * wz(ptsx-1,j) - 1070.d0 * wz(ptsx-2,j)
     &                 + 780.d0 * wz(ptsx-3,j) -  305.d0 * wz(ptsx-4,j)
     &                 +  50.d0 * wz(ptsx-5,j) ) / 225.d0
          th(ptsx,j) = ( 770.d0 * th(ptsx-1,j) - 1070.d0 * th(ptsx-2,j)
     &                 + 780.d0 * th(ptsx-3,j) -  305.d0 * th(ptsx-4,j)
     &                 +  50.d0 * th(ptsx-5,j) ) / 225.d0
              end if
            end do
            call loop(1d-7)

            write(*,*) my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax)
            if (mod(t, 5000) .eq. 0) call escreve(t)

            if (mod(t, 100) .eq. 0) then
              maxlocal_wz = maxval(dabs(wz1(:,:) - wz(:,:)))
              maxlocal_th = maxval(dabs(th1(:,:) - th(:,:)))
              write(*,*) my_rank, t, ux(ptsx,jmax), uy(ptsx,jmax),
     &                   maxlocal_wz, maxlocal_th
              call isstat(maxlocal_wz)
              call isstat(maxlocal_th)
              if ( (maxlocal_wz .lt. 1d-8) .and. (t .ge. stpp) ) then
                write(*,*) 'Stationary state.'
                call escreve_theta
                return
              end if
            end if

c            write(*,*) my_rank, t
          end do
          call escreve_theta
      end select

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init(dv1, wz1)

      ! initialize the program main variables
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include '../comm.multi'
      include '../comm.fs'
      include 'mpif.h'
      integer i, j
      character*20 nome
      real*8 m, dv1(ptsx,jmax), wz1(ptsx,jmax),
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax), xad,
     &       ue, afil(ptsx), bfil(ptsx), cfil(ptsx), ueptsx(ptsx)
      common/filt/ afil, bfil, cfil

      ! mount the lu matrix for the filter
      call lhs_tridf(afil, bfil, cfil)

      ! all the variables are set to zero
      do i = 1, ptsx
        do j = 1, jmax
          ux(i,j)  = 0.d0
          uy(i,j)  = 0.d0
          wz(i,j)  = 0.d0
          wz1(i,j) = 0.d0
          dv1(i,j) = 0.d0
        end do
      end do

      ! read the boundary layer profile
      open(1, file = '../pre_processing/base_fs.bin',
     &     form = 'unformatted')
      read(1) uxbt, uybt, wzbt
      close(unit = 1)

      ! give the values of the boundary layer profile for each node
      do j = 1, jmax
        do i = 1, ptsx
          ux(i,j) = uxbt(i+shift,j)
          uy(i,j) = uybt(i+shift,j)
          wz(i,j) = wzbt(i+shift,j)
        end do
      end do

      ! read the derivative and Poisson coefficients
      open(1, file = '../pre_processing/coefs.bin',
     &     form = 'unformatted')
      read(1) fp_fd_coef_e
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) sp_integ_coef
      read(1) cp_integ_coef
      read(1) pp_integ_coef
      read(1) lp_integ_coef
      close(unit = 1)

      ! mount the lhs for the derivative calculation
      call derivs_k

      ! read beta_fs from a file
      open(1, file = '../beta_fs.dist', form = 'formatted')
      read(1,*) beta_fs
      close(unit = 1)

      ! velocity at the edge of the boundary layer
      if (my_rank .eq. 0) then
        ue = ux(1,jmax)
      end if
      call MPI_BCAST(ue, 1, mpi_double_precision, 0, mpi_comm_world,
     &               ierr)

      ! define the boundary layer parameters
      do i = 1, ptsx
        xad        = x0 + dble(i+shift-1) * dx
        m          = beta_fs(i+shift) / (2.d0 - beta_fs(i+shift))
        ux(i,jmax) = ue*xad**m
        ueptsx(i)  = ux(i,jmax)
      end do
      call derparxue(duexmdx, ueptsx)

      ! calculate contants for multigrid method
      call create_ctes

      if (start .eq. 1) then
        write(nome,'(a,i0.2,a)')'based_',my_rank,'.bin'
        open (1, file = nome, form = 'unformatted')
        read(1) ux, uy, wz
        close(unit = 1)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init_theta(dv1, wz1, dvt1, th1)

      ! initialize the program main variables
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.coef'
      include '../comm.multi'
      include '../comm.fs'
      include 'mpif.h'
      integer i, j
      character*20 nome
      real*8 m, dv1(ptsx,jmax), wz1(ptsx,jmax),
     &       dvt1(ptsx,jmax), th1(ptsx,jmax), thbt(imax,jmax),
     &       uxbt(imax,jmax), uybt(imax,jmax), wzbt(imax,jmax), xad,
     &       ue, afil(ptsx), bfil(ptsx), cfil(ptsx), ueptsx(ptsx)

      common/filt/ afil, bfil, cfil
      ! mount the lu matrix for the filter
      call lhs_tridf(afil, bfil, cfil)

      ! all the variables are set to zero
      do i = 1, ptsx
        do j = 1, jmax
          ux(i,j)   = 0.d0
          uy(i,j)   = 0.d0
          wz(i,j)   = 0.d0
          th(i,j)   = 0.d0
          wz1(i,j)  = 0.d0
          dv1(i,j)  = 0.d0
          dvt1(i,j) = 0.d0
          th1(i,j)  = 0.d0
        end do
      end do

      ! reads the boundary layer profile
      open(1, file = '../pre_processing/base_fs.bin',
     &     form = 'unformatted')
      read(1) uxbt, uybt, wzbt, thbt
      close(unit = 1)

      ! give the values of the boundary layer profile for each node
      do j = 1, jmax
        do i = 1, ptsx
          ux(i,j) = uxbt(i+shift,j)
          uy(i,j) = uybt(i+shift,j)
          wz(i,j) = wzbt(i+shift,j)
          th(i,j) = thbt(i+shift,j)
        end do
      end do
 
      ! reads the derivative and Poisson coefficients
      open(1, file = '../pre_processing/coefs.bin',
     &     form = 'unformatted')
      read(1) fp_fd_coef_e
      read(1) sp_fd_coef
      read(1) cp_fd_coef
      read(1) pp_fd_coef
      read(1) lp_fd_coef
      read(1) fp_sd_coef
      read(1) sp_sd_coef
      read(1) cp_sd_coef
      read(1) pp_sd_coef
      read(1) lp_sd_coef
      read(1) sp_poi_coef
      read(1) cp_poi_coef
      read(1) pp_poi_coef
      read(1) lp_poi_coef
      read(1) w_at_w_coef
      read(1) dwydy_coef
      read(1) sp_integ_coef
      read(1) cp_integ_coef
      read(1) pp_integ_coef
      read(1) lp_integ_coef
      close(unit = 1)

      ! mount the lhs for the derivative calculation
      call derivs_k

      ! read beta_fs from a file
      open(1, file = '../beta_fs.dist', form = 'formatted')
      read(1,*) beta_fs
      close(unit = 1)

      ! velocity at the edge of the boundary layer
      if (my_rank .eq. 0) then
        ue = ux(1,jmax)
      end if
      call MPI_BCAST(ue, 1, mpi_double_precision, 0, mpi_comm_world,
     &               ierr)

      ! define the boundary layer parameters
      do i = 1, ptsx
        xad        = x0 + dble(i+shift-1) * dx
        m          = beta_fs(i+shift) / (2.d0 - beta_fs(i+shift))
        ux(i,jmax) = ue*xad**m
        ueptsx(i)  = ux(i,jmax)
      end do
      call derparxue(duexmdx, ueptsx)
      
      ! calculate contants for multigrid method
      call create_ctes

      if (start .eq. 1) then
        write(nome,'(a,i0.2,a)')'based_',my_rank,'.bin'
        open (1, file = nome, form = 'unformatted')
        read(1) ux, uy, wz, th
        close(unit = 1)
      end if


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine drv(dv)

      ! calculate the derivatives for RK method
      implicit none
      include '../par.for'
      include 'comm.var'
      include '../comm.par'
      integer i, j
      real*8 d2wzdx2(ptsx,jmax), d2wzdy2(ptsx,jmax),
     &       uwz(ptsx,jmax),     vwz(ptsx,jmax),
     &       duwzdx(ptsx,jmax),  dvwzdy(ptsx,jmax),
     &       uu(ptsx,jmax),      duudx(ptsx,jmax),    
     &       dv(ptsx,jmax), c(ptsx,jmax), erro

      do j = 1, jmax
        do i = 1, ptsx
          uwz(i,j) = ux(i,j) * wz(i,j)
          vwz(i,j) = uy(i,j) * wz(i,j)
          uu(i,j) = ux(i,j) * ux(i,j)
        end do
      end do

      ! x-direction derivatives
      call derparx(duwzdx, uwz)
      call derparx(duudx, uu)
      call derparxx(d2wzdx2, wz)

      ! y-direction derivatives
      call dery(dvwzdy, vwz)
      call deryy(d2wzdy2, wz)

      ! Ensuring outflow BC
      if (my_rank .eq. numproc) then
        d2wzdx2(ptsx,:) = 0.d0
      end if

      ! RHS of RK method
      do j = 2, jmax
        do i = 1, ptsx
          dv(i,j) = - duwzdx(i,j) - dvwzdy(i,j) 
     &              + Go2Re_base*duudx(i,j)
     &              + ( d2wzdx2(i,j) + fac_y * d2wzdy2(i,j) ) / Re
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine drv_theta(dv, dvt)

      ! calculate the derivatives for RK method
      implicit none
      include '../par.for'
      include 'comm.var'
      include '../comm.par'
      integer i, j
      real*8 d2wzdx2(ptsx,jmax), d2wzdy2(ptsx,jmax),
     &           uwz(ptsx,jmax),     vwz(ptsx,jmax),
     &        duwzdx(ptsx,jmax),  dvwzdy(ptsx,jmax),
     &            dv(ptsx,jmax),     uth(ptsx,jmax),
     &           vth(ptsx,jmax), d2thdx2(ptsx,jmax),
     &       d2thdy2(ptsx,jmax),     dvt(ptsx,jmax),
     &        duthdx(ptsx,jmax),  dvthdy(ptsx,jmax),
     &        uu(ptsx,jmax), duudx(ptsx,jmax)

      do j = 1, jmax
        do i = 1, ptsx
          uwz(i,j) = ux(i,j) * wz(i,j)
          vwz(i,j) = uy(i,j) * wz(i,j)
          uth(i,j) = ux(i,j) * th(i,j)
          vth(i,j) = uy(i,j) * th(i,j)
          uu(i,j)  = ux(i,j) * ux(i,j)
        end do
      end do

      ! first derivatives in x-direction
      call derparx(duwzdx, uwz)
      call derparx(duthdx, uth)
      call derparx(duudx, uu)

      ! first derivatives in y-direction
      call dery(dvwzdy, vwz)
      call dery(dvthdy, vth)

      ! second derivatives in x-direction
      call derparxx(d2wzdx2, wz)
      call derparxx(d2thdx2, th)

      ! second derivatives in y-direction
      call deryy(d2wzdy2, wz)
      call deryy(d2thdy2, th)

      ! Ensuring outflow BC
      if (my_rank .eq. numproc) then
        d2wzdx2(ptsx,:) = 0.d0
c        d2thdx2(ptsx,:) = 0.d0
      end if

      do j = 2, jmax
        do i = 1, ptsx
          dv(i,j)  = - duwzdx(i,j) - dvwzdy(i,j)
     &               + Go2Re_base * duudx(i,j)
     &               + ( d2wzdx2(i,j) + fac_y * d2wzdy2(i,j) ) / Re
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%% VER SE O FAC_Y ENTRA NESTA EQUAÇÃO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          dvt(i,j) = - duthdx(i,j) - dvthdy(i,j)
     &               + ( d2thdx2(i,j) + fac_y * d2thdy2(i,j) ) 
     &               / (Re * Pr)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine loop(erro)

      ! update the variables for each RK-step
      implicit none
      include '../par.for'
      real*8 erro

      if (erro .lt. 1d-4) call filter_trid
      call outuy
      call poi_uy(erro)
      call poi_ux
      call wz_wall

      return
      end
           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine create_ctes

      implicit none
      include '../par.for'
      include '../comm.multi'
      integer lvl, j
      real*8 aux
 
      ! multigrid spatial calculations
 
      ! dy0 at each multigrid level
      do lvl = 1 , msh
        if (stf .ne. 1.d0) then
          v_dy0(lvl) = dy * ( ( stf**(2**(lvl-1)) - 1.d0)
     &               / (stf - 1.d0) )
         else
          v_dy0(lvl) = dy * 2.d0**(lvl-1)
        end if
      end do

      do lvl = 1 , msh
        v_stf(lvl) = stf ** ( 2**(lvl-1) )
        v_dx2(lvl) = (dx * dble(2**(lvl-1)))**2
      end do

      ! dy at each space
      if (stf .ne. 1.d0) then
        do j = 1 , jmax - 1
          v_qdy(j)  = 1.d0 / (v_dy0(1) * v_stf(1)**(j-1))
        end do
       else
        do j = 1 , jmax - 1
          v_qdy(j)  = 1.d0 / (v_dy0(1))
        end do
      end if

      do lvl = 1 , msh
        if (stf .ne. 1.d0) then
          do j = 1 , (jmax - 1) / 2**(lvl-1)
            aux           = (v_dy0(lvl) * v_stf(lvl)**(j-1))
            v_qdy2(j,lvl) = 1.d0 / (aux**2)
          end do
         else
          do j = 1 , (jmax - 1) / 2**(lvl-1)
            aux           = (v_dy0(lvl))
            v_qdy2(j,lvl) = 1.d0 / (aux**2)
          end do
        end if
      end do

      v_ptsx(1) = ptsx
      v_ptsy(1) = jmax
      do lvl = 2 , msh
        v_ptsx(lvl) = (v_ptsx(lvl-1) + 1) / 2
        v_ptsy(lvl) = (v_ptsy(lvl-1) + 1) / 2
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine isstat(value_global)

      ! verifies if the current simulation has converged in each node
      ! and gives the answer to stop or continue
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer it
      integer status(MPI_status_size)
      real*8 value_global, value_local

      if (my_rank .gt. 0) then
        ! Send the error to the node 1
        call MPI_Send(value_global, 1, MPI_DOUBLE_PRECISION, 0, 51,
     &                MPI_COMM_WORLD, ierr)
        ! Receive answer if continues or not
        call MPI_Recv(value_global, 1, MPI_DOUBLE_PRECISION, 0, 71,
     &                MPI_COMM_WORLD, status, ierr)
       else
        do it = 1, numproc
          ! Receive the error to end or continue the program from 
          ! nodes 1 to numproc
          call MPI_Recv(value_local, 1, MPI_DOUBLE_PRECISION, it, 51,
     &                  MPI_COMM_WORLD, status, ierr)
          value_global = max(value_global,value_local)
c         value_global = value_global + value_local
        end do
        do it = 1, numproc
          call MPI_Send(value_global, 1, MPI_DOUBLE_PRECISION, it, 71,
     &                  MPI_COMM_WORLD,ierr)
        end do
        write(*,*) 'SC', value_global 
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve(t)

      ! write the results in binary form
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'

      integer t, c
      character*25 nm

      c = t / 1000

      write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
c     write(nm,'(a,i0.2,a,i0.3,a)')'based_',my_rank,'_',c,'.bin'
      open(1, file = nm, form = 'unformatted')
      write(1) ux, uy, wz
      close (unit = 1)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve2(t)

      ! write the results in binary form
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'

      integer t
      character*25 nm

      write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
      open(1, file = nm, form = 'unformatted')
      write(1) ux, uy, wz
      close (unit = 1)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine escreve_theta

      ! write the results in binary form
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      character*20 nm
 
      write(nm,'(a,i0.2,a)')'based_',my_rank,'.bin'
      open(1, file = nm, form = 'unformatted')
      write(1) ux, uy, wz, th
      close (unit = 1)

      return
      end 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                         end of main subroutines                              c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
