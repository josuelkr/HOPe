cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                          Poisson solver subroutine                           c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine poi_uy(erromax)

      ! parallel multigrid Poisson solver O(dx4,dy6) for the v-Poisson
      ! it makes v cicles until the defect is lower then the
      ! prescribed error (erromax)
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'comm.var'
      include '../comm.multi'
      logical continua
      integer i, j, itera, i_ini, i_fim, mesh, vcycles
      real*8 sum_err, erromax, fc(ptsx,jmax,msh), sc(ptsx,jmax,msh),
     &       fct(ptsx,jmax,msh)

      ! number of iterations in the SOR method
      itera = 2

      ! first point of the domain
      i_ini = 1
      if (my_rank .eq. 0) i_ini = 2

      ! last point of the domain
      i_fim = ptsx
      if (my_rank .eq. numproc) i_fim = ptsx - 1

      call function_source(fc, sc)
      continua = .true.
      do vcycles = 1, 5
        ! from finest to coarsest mesh
        do mesh = 1, msh-1
          call sor(fc, sc, itera, mesh)
          call defect(fc, sc, fct, -1.d0, mesh)
          ! stop condition
          if (mesh.eq.1) then
            call verifica_poi(fct, erromax, sum_err, continua)
            if (.not. continua) goto 22
          end if
          call restrict(fc, sc, fct, mesh)
          call defect(fc, sc, fct, 1.d0, mesh+1)
        end do
        ! coarsest mesh
        call sor(fc, sc, 100, msh)
        ! from coarsest to finest mesh
        do mesh = msh, 2, -1
          call correction(fc, fct, mesh)
          call intpol(fc, mesh)
          call sor(fc, sc, 1, mesh-1)
        end do                ! from coarsest to finest mesh >>>>>
      end do
   22 call boundary_exchange_end(fc)
      do j = 2, jmax
        do i = i_ini, i_fim
          uy(i,j) = fc(i,j,1)
        end do
      end do
      if (my_rank .eq. 0) then
        write(*,*) vcycles - 1, sum_err
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine verifica_poi(fct, errmax, sum_err, continua)

      ! verifies if the problem has converged in each node
      ! and gives the answer to stop or continue
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      logical continua
      integer it
      integer status(MPI_status_size)
      real*8 err, sum_err, errmax, fct(ptsx,jmax,msh)

      call get_max_value(fct, err)

      if (my_rank .gt. 0) then
        ! send the error to the node 1
        call MPI_Send(err,1, MPI_DOUBLE_PRECISION, 0, 51,
     &                MPI_COMM_WORLD, ierr)
        ! receive answer if continues or not
        call MPI_Recv(continua, 1, MPI_LOGICAL, 0, 71,
     &                MPI_COMM_WORLD, status, ierr)
       else
        sum_err = err
        do it = 1, numproc
          ! receive the error to end or continue the program from
          ! nodes 1 to numproc
          call MPI_Recv(err, 1, MPI_DOUBLE_PRECISION, it, 51,
     &                  MPI_COMM_WORLD, status, ierr)
          sum_err = max(sum_err,err)
        end do
        if (sum_err .lt. errmax) continua = .false.
        ! send answer if continues or not
        do it = 1, numproc
          call MPI_Send(continua, 1, MPI_LOGICAL, it, 71,
     &                  MPI_COMM_WORLD,ierr)
        end do
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine function_source(fc, sc)

      ! function and source terms of the finest grid Poisson equation
      implicit none
      include '../par.for'
      include 'comm.var'
      include '../comm.coef'
      integer i, j
      real*8 fc(ptsx,jmax,msh), sc(ptsx,jmax,msh),
     &       dwzdx(ptsx,jmax), dya
      common/dwdx/ dwzdx

      ! source term of the Poisson Equation
      do i = 1, ptsx
        ! weight near the boundaries
        sc(i,1,1) = 0.d0
        sc(i,2,1) = - sp_poi_coef(2,1) * dwzdx(i,2)
     &              - sp_poi_coef(3,1) * dwzdx(i,3)
        ! weight in the middle
        do j = 3, jmax - 2
          sc(i,j,1) = - cp_poi_coef(1,1) * dwzdx(i,j-1)
     &                - cp_poi_coef(2,1) * dwzdx(i,j)
     &                - cp_poi_coef(3,1) * dwzdx(i,j+1)
        end do
        ! weight near the boundaries
        sc(i,jmax-1,1) = - dwzdx(i,jmax-1)
        dya = dy * stf**(jmax-3)
        sc(i,jmax,1) = - lp_poi_coef(1,1) * dwzdx(i,jmax)
     &                 + lp_poi_coef(2,1) * duexmdx(i) / dya
      end do

      ! function of the Poisson equation
      do j = 1, jmax
        do i = 1, ptsx
          fc(i,j,1) = uy(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sor(fc, sc, itera, mesh)

      ! LSOR solver for each mesh of Poisson equation
      implicit none
      include '../par.for'
      include '../comm.par'
      include '../comm.multi'
      integer i, j, it, itera, mesh, indy(jmax), indyb(jmax)
      real*8 rf, ac(jmax,5), alc(jmax,5), acb(jmax,5), alcb(jmax,5),
     &       fc(ptsx,jmax,msh), sc(ptsx,jmax,msh), rhs(jmax)

      rf = 1.2d0
      if (itera .eq. 1) rf = 1.d0

      ! LHS calculation
      call coef_center(ac, mesh)
      call band5_poi(ac, v_ptsy(mesh), alc, indy)

      call coef_border(acb, mesh)
      call band5_poi(acb, v_ptsy(mesh), alcb, indyb)

      do it = 1, itera
        call boundary_exchange(fc, v_ptsx(mesh), v_ptsy(mesh), mesh)

        if (my_rank .eq. 0) then
          i = 2
          call rhs_border1(fc, sc, rhs, i, mesh)
          call banbk5_poi(acb, v_ptsy(mesh), alcb, indyb, rhs)
          do j = 2, v_ptsy(mesh)
            fc(i,j,mesh) = fc(i,j,mesh) + rf * (rhs(j) - fc(i,j,mesh))
          end do
        end if

        if (my_rank .eq. numproc) then
          i = v_ptsx(mesh) - 1
          call rhs_border2(fc, sc, rhs, i, mesh)
          call banbk5_poi(acb, v_ptsy(mesh), alcb, indyb, rhs)
          do j = 2, v_ptsy(mesh)
            fc(i,j,mesh) = fc(i,j,mesh) + rf * (rhs(j) - fc(i,j,mesh))
          end do
        end if

        do i = 3, v_ptsx(mesh) - 2
          call rhs_center(fc, sc, rhs, i, mesh)
          call banbk5_poi(ac, v_ptsy(mesh), alc, indy, rhs)
          do j = 2, v_ptsy(mesh)
            fc(i,j,mesh) = fc(i,j,mesh) + rf * (rhs(j) - fc(i,j,mesh))
          end do
        end do
      end do
    
      call boundary_exchange(fc, v_ptsx(mesh), v_ptsy(mesh), mesh)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coef_center(a, mesh)

      ! gives the LHS of the pentadiagonal matrix for the center
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer j, mesh
      real*8 a(jmax,5), qx

      qx = - 30.d0 / ( 12.d0 * v_dx2(mesh) * fac_y )

      a(1,1)   = 0.d0
      a(1,2)   = 0.d0
      a(1,3)   = 1.d0
      a(1,4)   = 0.d0
      a(1,5)   = 0.d0

      a(2,1)   = 0.d0
      a(2,2)   = sp_poi_coef(4,mesh) * v_qdy2(1,mesh)
      a(2,3)   = sp_poi_coef(5,mesh) * v_qdy2(1,mesh)
     &         + sp_poi_coef(2,mesh) * qx
      a(2,4)   = sp_poi_coef(6,mesh) * v_qdy2(1,mesh)
     &         + sp_poi_coef(3,mesh) * qx
      a(2,5)   = sp_poi_coef(7,mesh) * v_qdy2(1,mesh)

      do j = 3, v_ptsy(mesh) - 2
        a(j,1) = cp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
        a(j,2) = cp_poi_coef(5,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(1,mesh) * qx
        a(j,3) = cp_poi_coef(6,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(2,mesh) * qx
        a(j,4) = cp_poi_coef(7,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(3,mesh) * qx
        a(j,5) = cp_poi_coef(8,mesh) * v_qdy2(j-2,mesh)
      end do

      j = v_ptsy(mesh)
      a(j-1,1) = 0.d0
      a(j-1,2) = pp_poi_coef(2,mesh) * v_qdy2(j-2,mesh)
      a(j-1,3) = pp_poi_coef(3,mesh) * v_qdy2(j-2,mesh) + qx
      a(j-1,4) = pp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
      a(j-1,5) = 0.d0
 
      a(j,1)   = lp_poi_coef(5,mesh) * v_qdy2(j-2,mesh)
      a(j,2)   = lp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
      a(j,3)   = lp_poi_coef(3,mesh) * v_qdy2(j-2,mesh) + qx
      a(j,4)   = 0.d0
      a(j,5)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coef_border(a, mesh)

      ! gives the LHS of the pentadiagonal matrix for the borders
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer j, mesh
      real*8 a(jmax,5), qx

      qx = - 15.d0  / ( 12.d0 * v_dx2(mesh) * fac_y )

      a(1,1)   = 0.d0
      a(1,2)   = 0.d0
      a(1,3)   = 1.d0
      a(1,4)   = 0.d0
      a(1,5)   = 0.d0

      a(2,1)   = 0.d0
      a(2,2)   = sp_poi_coef(4,mesh) * v_qdy2(1,mesh)
      a(2,3)   = sp_poi_coef(5,mesh) * v_qdy2(1,mesh)
     &         + sp_poi_coef(2,mesh) * qx
      a(2,4)   = sp_poi_coef(6,mesh) * v_qdy2(1,mesh)
     &         + sp_poi_coef(3,mesh) * qx
      a(2,5)   = sp_poi_coef(7,mesh) * v_qdy2(1,mesh)

      do j = 3, v_ptsy(mesh) - 2
       a(j,1)  = cp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
       a(j,2)  = cp_poi_coef(5,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(1,mesh) * qx
       a(j,3)  = cp_poi_coef(6,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(2,mesh) * qx
       a(j,4)  = cp_poi_coef(7,mesh) * v_qdy2(j-2,mesh)
     &         + cp_poi_coef(3,mesh) * qx
       a(j,5)  = cp_poi_coef(8,mesh) * v_qdy2(j-2,mesh)
      end do

      j = v_ptsy(mesh)
      a(j-1,1) = 0.d0
      a(j-1,2) = pp_poi_coef(2,mesh) * v_qdy2(j-2,mesh)
      a(j-1,3) = pp_poi_coef(3,mesh) * v_qdy2(j-2,mesh) + qx
      a(j-1,4) = pp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
      a(j-1,5) = 0.d0
 
      a(j,1)   = lp_poi_coef(5,mesh) * v_qdy2(j-2,mesh)
      a(j,2)   = lp_poi_coef(4,mesh) * v_qdy2(j-2,mesh)
      a(j,3)   = lp_poi_coef(3,mesh) * v_qdy2(j-2,mesh) + qx
      a(j,4)   = 0.d0
      a(j,5)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhs_border1(fc, sc, rhs, i, mesh)

      ! RHS for the first derivative calculation in y direction
      ! used just by my_rank=0
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j, mesh
      real*8 qx, rhs(jmax), fc(ptsx,jmax,msh), sc(ptsx,jmax,msh)

      qx = 1.d0 / ( 12.d0 *  v_dx2(mesh) * fac_y )

      rhs(1) = fc(i,1,mesh)

      j = 2
      rhs(j) = sc(i,j,mesh) - qx
     &       * ( sp_poi_coef(3,mesh) * ( 10.d0 * fc(i-1,j+1,mesh)
     &                                 -  4.d0 * fc(i+1,j+1,mesh)
     &                                 + 14.d0 * fc(i+2,j+1,mesh)
     &                                 -  6.d0 * fc(i+3,j+1,mesh)
     &                                 +         fc(i+4,j+1,mesh) )
     &         + sp_poi_coef(2,mesh) * ( 10.d0 * fc(i-1,j,mesh)
     &                                 -  4.d0 * fc(i+1,j,mesh)
     &                                 + 14.d0 * fc(i+2,j,mesh)
     &                                 -  6.d0 * fc(i+3,j,mesh)
     &                                 +         fc(i+4,j,mesh) ) )

      do j = 3, v_ptsy(mesh) - 2
        rhs(j) = sc(i,j,mesh) - qx
     &         * ( cp_poi_coef(1,mesh) * ( 10.d0 * fc(i-1,j-1,mesh)
     &                                   -  4.d0 * fc(i+1,j-1,mesh)
     &                                   + 14.d0 * fc(i+2,j-1,mesh)
     &                                   -  6.d0 * fc(i+3,j-1,mesh)
     &                                   +         fc(i+4,j-1,mesh) )
     &           + cp_poi_coef(2,mesh) * ( 10.d0 * fc(i-1,j,mesh)
     &                                   -  4.d0 * fc(i+1,j,mesh)
     &                                   + 14.d0 * fc(i+2,j,mesh)
     &                                   -  6.d0 * fc(i+3,j,mesh)
     &                                   +         fc(i+4,j,mesh) )
     &           + cp_poi_coef(3,mesh) * ( 10.d0 * fc(i-1,j+1,mesh)
     &                                   -  4.d0 * fc(i+1,j+1,mesh)
     &                                   + 14.d0 * fc(i+2,j+1,mesh)
     &                                   -  6.d0 * fc(i+3,j+1,mesh)
     &                                   +         fc(i+4,j+1,mesh) ) )
      end do

      j = v_ptsy(mesh) - 1
      rhs(j) = sc(i,j,mesh) - qx * ( 10.d0 * fc(i-1,j,mesh)
     &                             -  4.d0 * fc(i+1,j,mesh)
     &                             + 14.d0 * fc(i+2,j,mesh)
     &                             -  6.d0 * fc(i+3,j,mesh)
     &                             +         fc(i+4,j,mesh) )

      j = j + 1
      rhs(j) = sc(i,j,mesh) - qx * ( 10.d0 * fc(i-1,j,mesh)
     &                             -  4.d0 * fc(i+1,j,mesh)
     &                             + 14.d0 * fc(i+2,j,mesh)
     &                             -  6.d0 * fc(i+3,j,mesh)
     &                             +         fc(i+4,j,mesh) )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhs_border2(fc, sc, rhs, i, mesh)

      ! RHS for the first derivative calculation in y direction
      ! used just by my_rank=numproc
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j, mesh
      real*8 qx, rhs(jmax), fc(ptsx,jmax,msh), sc(ptsx,jmax,msh)

      qx = 1.d0 / ( 12.d0 * v_dx2(mesh) * fac_y )

      rhs(1) = fc(i,1,mesh)

      j = 2
      rhs(j) = sc(i,j,mesh) - qx
     &       * ( sp_poi_coef(3,mesh) * ( 10.d0 * fc(i+1,j+1,mesh)
     &                                 -  4.d0 * fc(i-1,j+1,mesh)
     &                                 + 14.d0 * fc(i-2,j+1,mesh)
     &                                 -  6.d0 * fc(i-3,j+1,mesh)
     &                                 +         fc(i-4,j+1,mesh) )
     &         + sp_poi_coef(2,mesh) * ( 10.d0 * fc(i+1,j,mesh)
     &                                 -  4.d0 * fc(i-1,j,mesh)
     &                                 + 14.d0 * fc(i-2,j,mesh)
     &                                 -  6.d0 * fc(i-3,j,mesh)
     &                                 +         fc(i-4,j,mesh) ) )

      do j = 3, v_ptsy(mesh) - 2
        rhs(j) = sc(i,j,mesh) - qx
     &         * ( cp_poi_coef(1,mesh) * ( 10.d0 * fc(i+1,j-1,mesh)
     &                                   -  4.d0 * fc(i-1,j-1,mesh)
     &                                   + 14.d0 * fc(i-2,j-1,mesh)
     &                                   -  6.d0 * fc(i-3,j-1,mesh)
     &                                   +         fc(i-4,j-1,mesh) )
     &           + cp_poi_coef(2,mesh) * ( 10.d0 * fc(i+1,j,mesh)
     &                                   -  4.d0 * fc(i-1,j,mesh)
     &                                   + 14.d0 * fc(i-2,j,mesh)
     &                                   -  6.d0 * fc(i-3,j,mesh)
     &                                   +         fc(i-4,j,mesh) )
     &           + cp_poi_coef(3,mesh) * ( 10.d0 * fc(i+1,j+1,mesh)
     &                                   -  4.d0 * fc(i-1,j+1,mesh)
     &                                   + 14.d0 * fc(i-2,j+1,mesh)
     &                                   -  6.d0 * fc(i-3,j+1,mesh)
     &                                   +         fc(i-4,j+1,mesh) ) )
      end do

      j = v_ptsy(mesh) - 1
      rhs(j) = sc(i,j,mesh) - qx * ( 10.d0 * fc(i+1,j,mesh)
     &                             -  4.d0 * fc(i-1,j,mesh)
     &                             + 14.d0 * fc(i-2,j,mesh)
     &                             -  6.d0 * fc(i-3,j,mesh)
     &                             +         fc(i-4,j,mesh) )

      j = j + 1
      rhs(j) = sc(i,j,mesh) - qx * ( 10.d0 * fc(i+1,j,mesh)
     &                             -  4.d0 * fc(i-1,j,mesh)
     &                             + 14.d0 * fc(i-2,j,mesh)
     &                             -  6.d0 * fc(i-3,j,mesh)
     &                             +         fc(i-4,j,mesh) )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhs_center(fc, sc, rhs, i, mesh)

      ! RHS for the first derivative calculation in y direction
      ! used for all nodes
      implicit none
      include '../par.for'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j, mesh
      real*8 qx, rhs(jmax), fc(ptsx,jmax,msh), sc(ptsx,jmax,msh)
      
      qx = 1.d0 / ( 12.d0 * v_dx2(mesh) * fac_y )

      rhs(1) = fc(i,1,mesh)

      j = 2
      rhs(j) = sc(i,j,mesh) - qx
     &       * ( sp_poi_coef(3,mesh) * ( -  1.d0 * fc(i-2,j+1,mesh)
     &                                   + 16.d0 * fc(i-1,j+1,mesh)
     &                                   + 16.d0 * fc(i+1,j+1,mesh)
     &                                   -  1.d0 * fc(i+2,j+1,mesh) )
     &         + sp_poi_coef(2,mesh) * ( -  1.d0 * fc(i-2,j,mesh)
     &                                   + 16.d0 * fc(i-1,j,mesh)
     &                                   + 16.d0 * fc(i+1,j,mesh)
     &                                   -  1.d0 * fc(i+2,j,mesh) ) )

      do j = 3, v_ptsy(mesh) - 2
        rhs(j) = sc(i,j,mesh) - qx
     &         * ( cp_poi_coef(1,mesh)
     &         * ( -  1.d0 * fc(i-2,j-1,mesh)
     &             + 16.d0 * fc(i-1,j-1,mesh)
     &             + 16.d0 * fc(i+1,j-1,mesh)
     &             -  1.d0 * fc(i+2,j-1,mesh) )
     &           + cp_poi_coef(2,mesh)
     &         * ( -  1.d0 * fc(i-2,j,mesh)
     &             + 16.d0 * fc(i-1,j,mesh)
     &             + 16.d0 * fc(i+1,j,mesh)
     &             -  1.d0 * fc(i+2,j,mesh) )
     &           + cp_poi_coef(3,mesh)
     &         * ( -  1.d0 * fc(i-2,j+1,mesh)
     &             + 16.d0 * fc(i-1,j+1,mesh)
     &             + 16.d0 * fc(i+1,j+1,mesh)
     &             -  1.d0 * fc(i+2,j+1,mesh) ) )
      end do

      j = v_ptsy(mesh) - 1
      rhs(j) = sc(i,j,mesh) - qx
     &       * (       - ( fc(i+2,j,mesh) + fc(i-2,j,mesh) )
     &         + 16.d0 * ( fc(i+1,j,mesh) + fc(i-1,j,mesh) ) )

      j = j + 1
      rhs(j) = sc(i,j,mesh) - qx
     &       * (       - ( fc(i+2,j,mesh) + fc(i-2,j,mesh) )
     &         + 16.d0 * ( fc(i+1,j,mesh) + fc(i-1,j,mesh) ) )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine defect(fc, sc, fct, si, mesh)

      ! calculates the residue or the new source term of the
      ! Poisson equation depending on the si
      implicit none
      include '../par.for'
      include '../comm.par'
      include '../comm.coef'
      include '../comm.multi'
      integer i, j, mesh
      real*8  si, b, c, fc(ptsx,jmax,msh), sc(ptsx,jmax,msh),
     &        fct(ptsx,jmax,msh)

      if (my_rank.eq.0) then
        ! point ( i = 2, j = 2 )
        i = 2
        j = 2
        b = ( sp_poi_coef(3,mesh)
     &    * ( 10.d0 * fc(i-1,j+1,mesh) - 15.d0 * fc(i,j+1,mesh)
     &      -  4.d0 * fc(i+1,j+1,mesh) + 14.d0 * fc(i+2,j+1,mesh)
     &      -  6.d0 * fc(i+3,j+1,mesh) +         fc(i+4,j+1,mesh) )
     &      + sp_poi_coef(2,mesh)
     &    * ( 10.d0 * fc(i-1,j,mesh)   - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i+1,j,mesh)   + 14.d0 * fc(i+2,j,mesh)
     &      -  6.d0 * fc(i+3,j,mesh)   +         fc(i+4,j,mesh) ) )
        c = sp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &    + sp_poi_coef(5,mesh) * fc(i,j,mesh)
     &    + sp_poi_coef(6,mesh) * fc(i,j+1,mesh)
     &    + sp_poi_coef(7,mesh) * fc(i,j+2,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(1,mesh) )
        ! middle of the second column ( i = 2, 3>=j>=nptsy-2 )
        do j = 3, v_ptsy(mesh) - 2
          b = ( cp_poi_coef(1,mesh)
     &      * ( 10.d0 * fc(i-1,j-1,mesh) - 15.d0 * fc(i,j-1,mesh)
     &        -  4.d0 * fc(i+1,j-1,mesh) + 14.d0 * fc(i+2,j-1,mesh)
     &        -  6.d0 * fc(i+3,j-1,mesh) +         fc(i+4,j-1,mesh) )
     &        + cp_poi_coef(2,mesh)
     &      * ( 10.d0 * fc(i-1,j,mesh)   - 15.d0 * fc(i,j,mesh)
     &        -  4.d0 * fc(i+1,j,mesh)   + 14.d0 * fc(i+2,j,mesh)
     &        -  6.d0 * fc(i+3,j,mesh)   +         fc(i+4,j,mesh) )
     &        + cp_poi_coef(3,mesh)
     &      * ( 10.d0 * fc(i-1,j+1,mesh) - 15.d0 * fc(i,j+1,mesh)
     &        -  4.d0 * fc(i+1,j+1,mesh) + 14.d0 * fc(i+2,j+1,mesh)
     &        -  6.d0 * fc(i+3,j+1,mesh) +         fc(i+4,j+1,mesh) ) )
          c = ( cp_poi_coef(4,mesh) * fc(i,j-2,mesh)
     &        + cp_poi_coef(5,mesh) * fc(i,j-1,mesh)
     &        + cp_poi_coef(6,mesh) * fc(i,j,mesh)
     &        + cp_poi_coef(7,mesh) * fc(i,j+1,mesh)
     &        + cp_poi_coef(8,mesh) * fc(i,j+2,mesh) )
          fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &               * fac_y) + c * v_qdy2(j-2,mesh) )
        end do
        ! point ( i = 2, j = nptsy - 1 )
        j = v_ptsy(mesh) - 1
        b = ( 10.d0 * fc(i-1,j,mesh) - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i+1,j,mesh) + 14.d0 * fc(i+2,j,mesh)
     &      -  6.d0 * fc(i+3,j,mesh) +         fc(i+4,j,mesh) )
        c = pp_poi_coef(4,mesh) * fc(i,j+1,mesh)
     &    + pp_poi_coef(3,mesh) * fc(i,j,mesh)
     &    + pp_poi_coef(2,mesh) * fc(i,j-1,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-1,mesh) )
        ! point ( i = 2, j = nptsy )
        j = v_ptsy(mesh)
        b = ( 10.d0 * fc(i-1,j,mesh) - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i+1,j,mesh) + 14.d0 * fc(i+2,j,mesh)
     &      -  6.d0 * fc(i+3,j,mesh) +         fc(i+4,j,mesh) )
        c = lp_poi_coef(3,mesh) * fc(i,j,mesh)
     &    + lp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &    + lp_poi_coef(5,mesh) * fc(i,j-2,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-2,mesh) )
      end if

      if (my_rank .eq. numproc) then
        ! point ( i = nptsx - 1, j = 2 )
        i = v_ptsx(mesh) - 1
        j = 2
        b = sp_poi_coef(3,mesh)
     &    * ( 10.d0 * fc(i+1,j+1,mesh) - 15.d0 * fc(i,j+1,mesh)
     &      -  4.d0 * fc(i-1,j+1,mesh) + 14.d0 * fc(i-2,j+1,mesh)
     &      -  6.d0 * fc(i-3,j+1,mesh) +         fc(i-4,j+1,mesh) )
     &    + sp_poi_coef(2,mesh)
     &    * ( 10.d0 * fc(i+1,j,mesh)   - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i-1,j,mesh)   + 14.d0 * fc(i-2,j,mesh)
     &      -  6.d0 * fc(i-3,j,mesh)   +         fc(i-4,j,mesh) )
        c = sp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &    + sp_poi_coef(5,mesh) * fc(i,j,mesh)
     &    + sp_poi_coef(6,mesh) * fc(i,j+1,mesh)
     &    + sp_poi_coef(7,mesh) * fc(i,j+2,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(1,mesh) )
        ! middle of the last but one column ( i = nptsx - 1, 3>=j>=nptsy-2 )
        do j = 3, v_ptsy(mesh) - 2
          b = cp_poi_coef(1,mesh)
     &      * ( 10.d0 * fc(i+1,j-1,mesh) - 15.d0 * fc(i,j-1,mesh)
     &        -  4.d0 * fc(i-1,j-1,mesh) + 14.d0 * fc(i-2,j-1,mesh)
     &        -  6.d0 * fc(i-3,j-1,mesh) +         fc(i-4,j-1,mesh) )
     &      + cp_poi_coef(2,mesh)
     &      * ( 10.d0 * fc(i+1,j,mesh)   - 15.d0 * fc(i,j,mesh)
     &        -  4.d0 * fc(i-1,j,mesh)   + 14.d0 * fc(i-2,j,mesh)
     &        -  6.d0 * fc(i-3,j,mesh)   +         fc(i-4,j,mesh) )
     &      + cp_poi_coef(3,mesh)
     &      * ( 10.d0 * fc(i+1,j+1,mesh) - 15.d0 * fc(i,j+1,mesh)
     &        -  4.d0 * fc(i-1,j+1,mesh) + 14.d0 * fc(i-2,j+1,mesh)
     &        -  6.d0 * fc(i-3,j+1,mesh) +         fc(i-4,j+1,mesh) )
          c = cp_poi_coef(4,mesh) * fc(i,j-2,mesh)
     &      + cp_poi_coef(5,mesh) * fc(i,j-1,mesh)
     &      + cp_poi_coef(6,mesh) * fc(i,j,mesh)
     &      + cp_poi_coef(7,mesh) * fc(i,j+1,mesh)
     &      + cp_poi_coef(8,mesh) * fc(i,j+2,mesh)
          fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &               * fac_y) + c * v_qdy2(j-2,mesh) )
        end do
        ! point ( i = nptsx - 1, j = nptsy - 1 )
        j = v_ptsy(mesh) - 1
        b = ( 10.d0 * fc(i+1,j,mesh) - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i-1,j,mesh) + 14.d0 * fc(i-2,j,mesh)
     &      -  6.d0 * fc(i-3,j,mesh) +         fc(i-4,j,mesh) )
        c = pp_poi_coef(4,mesh) * fc(i,j+1,mesh)
     &    + pp_poi_coef(3,mesh) * fc(i,j,mesh)
     &    + pp_poi_coef(2,mesh) * fc(i,j-1,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-1,mesh) )
        ! point ( i = nptsx - 1, j = nptsy )
        j = v_ptsy(mesh)
        b = ( 10.d0 * fc(i+1,j,mesh) - 15.d0 * fc(i,j,mesh)
     &      -  4.d0 * fc(i-1,j,mesh) + 14.d0 * fc(i-2,j,mesh)
     &      -  6.d0 * fc(i-3,j,mesh) +         fc(i-4,j,mesh) )
        c = ( lp_poi_coef(3,mesh) * fc(i,j,mesh)
     &      + lp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &      + lp_poi_coef(5,mesh) * fc(i,j-2,mesh))
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-2,mesh) )
      end if

      ! middle region en x
      do i = 3, v_ptsx(mesh) - 2
        ! second line ( j = 2, 3>=i>=nptsx-2 )
        j = 2
        b = ( sp_poi_coef(3,mesh) * ( -  1.d0 * fc(i-2,j+1,mesh)
     &                                + 16.d0 * fc(i-1,j+1,mesh)
     &                                - 30.d0 * fc(i,j+1,mesh)
     &                                + 16.d0 * fc(i+1,j+1,mesh)
     &                                -  1.d0 * fc(i+2,j+1,mesh) )
     &      + sp_poi_coef(2,mesh) * ( -  1.d0 * fc(i-2,j,mesh)
     &                                + 16.d0 * fc(i-1,j,mesh)
     &                                - 30.d0 * fc(i,j,mesh)
     &                                + 16.d0 * fc(i+1,j,mesh)
     &                                -  1.d0 * fc(i+2,j,mesh) ) )
        c = sp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &    + sp_poi_coef(5,mesh) * fc(i,j,mesh)
     &    + sp_poi_coef(6,mesh) * fc(i,j+1,mesh)
     &    + sp_poi_coef(7,mesh) * fc(i,j+2,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(1,mesh) )
        ! middle region
        do j = 3, v_ptsy(mesh) - 2
          b = ( cp_poi_coef(1,mesh) * ( -  1.d0 * fc(i-2,j-1,mesh)
     &                                  + 16.d0 * fc(i-1,j-1,mesh)
     &                                  - 30.d0 * fc(i,j-1,mesh)
     &                                  + 16.d0 * fc(i+1,j-1,mesh)
     &                                  -  1.d0 * fc(i+2,j-1,mesh) )
     &        + cp_poi_coef(2,mesh) * ( -  1.d0 * fc(i-2,j,mesh)
     &                                  + 16.d0 * fc(i-1,j,mesh)
     &                                  - 30.d0 * fc(i,j,mesh)
     &                                  + 16.d0 * fc(i+1,j,mesh)
     &                                  -  1.d0 * fc(i+2,j,mesh) )
     &        + cp_poi_coef(3,mesh) * ( -  1.d0 * fc(i-2,j+1,mesh)
     &                                  + 16.d0 * fc(i-1,j+1,mesh)
     &                                  - 30.d0 * fc(i,j+1,mesh)
     &                                  + 16.d0 * fc(i+1,j+1,mesh)
     &                                  -  1.d0 * fc(i+2,j+1,mesh) ) )
          c = cp_poi_coef(4,mesh) * fc(i,j-2,mesh)
     &      + cp_poi_coef(5,mesh) * fc(i,j-1,mesh)
     &      + cp_poi_coef(6,mesh) * fc(i,j,mesh)
     &      + cp_poi_coef(7,mesh) * fc(i,j+1,mesh)
     &      + cp_poi_coef(8,mesh) * fc(i,j+2,mesh)
          fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &               * fac_y) + c * v_qdy2(j-2,mesh) )
        end do
        ! last but one line ( j = nptsy-1, 3>=i>=nptsx-2 )
        j = v_ptsy(mesh) - 1
        b = ( -  1.d0 * ( fc(i-2,j,mesh) + fc(i+2,j,mesh) )
     &        + 16.d0 * ( fc(i-1,j,mesh) + fc(i+1,j,mesh) )
     &        - 30.d0 *   fc(i,j,mesh) )
        c = pp_poi_coef(4,mesh) * fc(i,j+1,mesh)
     &    + pp_poi_coef(3,mesh) * fc(i,j,mesh)
     &    + pp_poi_coef(2,mesh) * fc(i,j-1,mesh)
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-1,mesh) )
        ! last line ( j = nptsy, 3>=i>=nptsx-2 )
        j = v_ptsy(mesh)
        b = ( -  1.d0 * ( fc(i-2,j,mesh) + fc(i+2,j,mesh) )
     &        + 16.d0 * ( fc(i-1,j,mesh) + fc(i+1,j,mesh) )
     &        - 30.d0 *   fc(i,j,mesh) )
        c = ( lp_poi_coef(3,mesh) * fc(i,j,mesh)
     &      + lp_poi_coef(4,mesh) * fc(i,j-1,mesh)
     &      + lp_poi_coef(5,mesh) * fc(i,j-2,mesh) )
        fct(i,j,1) = sc(i,j,mesh) + si * ( b / (12.d0 * v_dx2(mesh)
     &             * fac_y) + c * v_qdy2(j-2,mesh) )
      end do

      if (si .lt. 0) return

      do j = 2, v_ptsy(mesh)
        do i = 2, v_ptsx(mesh) - 1
          sc(i,j,mesh) = fct(i,j,1)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine restrict(fc, sc, fct, mesh)

      ! restricts the result from a fine mesh to a coarser one
      implicit none
      include '../par.for'
      include '../comm.par'
      include '../comm.multi'
      integer i, j, fi, fj, mesh
      real*8 fc(ptsx,jmax,msh), sc(ptsx,jmax,msh), fct(ptsx,jmax,msh)

      ! restriction in the middle
      fj = 3
      do j = 2, v_ptsy(mesh+1) - 1
        fi = 3
        do i = 2, v_ptsx(mesh+1) - 1
          fc(i,j,mesh+1)  = fc(fi,fj,mesh)
          fct(i,j,mesh+1) = fc(i,j,mesh+1)
          sc(i,j,mesh+1)  = 6.25d-2 * ( 2.d0 * ( fct(fi-1,fj,1)
     &                                         + fct(fi+1,fj,1)
     &                                         + fct(fi,fj-1,1)
     &                                         + fct(fi,fj+1,1) )
     &                                + 4.d0 *   fct(fi,fj,1)
     &                                         + fct(fi+1,fj+1,1)
     &                                         + fct(fi+1,fj-1,1)
     &                                         + fct(fi-1,fj+1,1)
     &                                         + fct(fi-1,fj-1,1) )
          fi = fi + 2
        end do
        fj = fj + 2
      end do

      ! restriction at upper and lower boundaries
      fi = 3
      do i = 2, v_ptsx(mesh+1) - 1
        j               = 1
        fc(i,j,mesh+1)  = fc(fi,j,mesh)
        fct(i,j,mesh+1) = fc(i,j,mesh+1)
        j               = v_ptsy(mesh+1)
        fj              = 2 * j - 1
        fc(i,j,mesh+1)  = fc(fi,fj,mesh)
        fct(i,j,mesh+1) = fc(i,j,mesh+1)
        sc(i,j,mesh+1)  = 0.25d0 * ( fct(fi-1,fj,1)
     &                  + 2.d0   *   fct(fi,fj,1)
     &                  +            fct(fi+1,fj,1) )
        fi              = fi + 2
      end do

      ! restriction in the first column
      if (my_rank .eq. 0) then
        i               = 1
        fc(i,1,mesh+1)  = fc(i,1,mesh)
        fct(i,1,mesh+1) = fc(i,1,mesh+1)
        fj              = 3
        do j = 2, v_ptsy(mesh+1) - 1
          fc(i,j,mesh+1)  = fc(i,fj,mesh)
          fct(i,j,mesh+1) = fc(i,j,mesh+1)
          fj              = fj + 2
        end do
        fc(i,v_ptsy(mesh+1),mesh+1)  = fc(i,v_ptsy(mesh),mesh)
        fct(i,v_ptsy(mesh+1),mesh+1) = fc(i,v_ptsy(mesh+1),mesh+1)
      end if

      ! restriction in the last column
      if (my_rank .eq. numproc) then
        i               = v_ptsx(mesh+1)
        fi              = 2 * i - 1
        fc(i,1,mesh+1)  = fc(fi,1,mesh)
        fct(i,1,mesh+1) = fc(i,1,mesh+1)
        fj              = 3
        do j = 2, v_ptsy(mesh+1) - 1
          fc(i,j,mesh+1)  = fc(fi,fj,mesh)
          fct(i,j,mesh+1) = fc(i,j,mesh+1)
          fj              = fj + 2
        end do
        fc(i,v_ptsy(mesh+1),mesh+1)  = fc(fi,v_ptsy(mesh),mesh)
        fct(i,v_ptsy(mesh+1),mesh+1) = fc(i,v_ptsy(mesh+1),mesh+1)
      end if

      call boundary_exchange(fc, v_ptsx(mesh+1), v_ptsy(mesh+1), mesh+1)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine intpol(fc, mesh)

      ! interpol the results from a coarse mesh to a fine one
      implicit none
      include '../par.for'
      include '../comm.par'
      include '../comm.multi'
      integer i, j, fi, fj, mesh
      real*8 fc(ptsx,jmax,msh), temp(ptsx,jmax)

      ! copying fcc to temp
      fj = 1
      do j = 1, v_ptsy(mesh)
        fi = 1
        do i = 1, v_ptsx(mesh)
          temp(fi,fj) = fc(i,j,mesh)
          fi = fi + 2
        end do
        fj = fj + 2
      end do

      ! Interpolation in x
      do j = 1, v_ptsy(mesh-1), 2
        do i = 2, v_ptsx(mesh-1) - 1, 2
          temp(i,j) = 0.5d0 * ( temp(i-1,j) + temp(i+1,j) )
        end do
      end do

      ! Interpolation in y
      do j = 2, v_ptsy(mesh-1) - 1, 2
        do i = 1, v_ptsx(mesh-1)
          temp(i,j) = 0.5d0 * ( temp(i,j+1) + temp(i,j-1) )
        end do
      end do

      ! adition to new function
      do j = 1, v_ptsy(mesh-1)
        do i = 2, v_ptsx(mesh-1) - 1
          fc(i,j,mesh-1) = fc(i,j,mesh-1) + temp(i,j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine correction(fc, fct, mesh)

      ! applies the correction when going from a coarse to fine mesh
      implicit none
      include '../par.for'
      include '../comm.par'
      include '../comm.multi'
      integer i, j, mesh
      real*8 fc(ptsx,jmax,msh), fct(ptsx,jmax,msh)

      do j = 1, v_ptsy(mesh)
        do i = 2, v_ptsx(mesh) - 1
          fc(i,j,mesh) = fc(i,j,mesh) - fct(i,j,mesh)
        end do
      end do

      call boundary_exchange(fc, v_ptsx(mesh), v_ptsy(mesh), mesh)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine get_max_value(fct, err_pc)

      implicit none
      include '../par.for'
      include 'comm.var'
      include '../comm.par'
      integer i, j
      real*8 err_pc, fct(ptsx,jmax,msh)

      err_pc = 0.d0

      do j = 2, jmax
        do i = 2, ptsx - 1
          err_pc = max( err_pc, abs(fct(i,j,1)) )
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange(var, nptsx, nptsy, mesh)

      ! exchange values of the boundaries
      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer status(MPI_Status_size)
      integer i, j, nptsy, nptsx, interm, varint, mesh
      real*8 aux(meshpx*nptsy), var(ptsx,jmax,msh)

      ! variable used to calculate the number of columns needed 
      ! to go forward or backward (interm)
      varint = ( ptsx - 1 ) / ( nptsx - 1)
      interm = inter / varint

      if (my_rank .lt. numproc) then
        ! Sending the near right boundary columns to node + 1
        do i = 1, meshpx
          do j = 1, nptsy
            aux(j+(i-1)*nptsy) = var(i + nptsx - interm - 1,j,mesh)
          end do
        end do
        call MPI_Send(aux, meshpx*nptsy, mpi_double_precision,
     &                my_rank + 1, 10, MPI_COMM_WORLD, ierr)
        ! Receiving the new right boundary columns from node + 1
        call MPI_Recv(aux, meshpx*nptsy, mpi_double_precision,
     &                my_rank + 1, 20, MPI_COMM_WORLD, status, ierr)
        do i = 1, meshpx
          do j = 1, nptsy
            var(i + nptsx - meshpx,j,mesh) = aux(j+(i-1)*nptsy)
          end do
        end do
      end if

      if (my_rank .gt. 0) then
        ! Receiving the new left boundary columns from node - 1
        call MPI_Recv(aux, meshpx*nptsy, mpi_double_precision,
     &                my_rank - 1, 10, MPI_COMM_WORLD, status, ierr)
        do i = 1, meshpx
          do j = 1, nptsy
            var(i,j,mesh) = aux(j+(i-1)*nptsy)
          end do
        end do
        ! Sending the new right boundary columns to node - 1
        do i = 1, meshpx
          do j = 1, nptsy
            aux(j+(i-1)*nptsy) = var(i + interm - meshpx + 1,j,mesh)
          end do
        end do
        call MPI_Send(aux, meshpx*nptsy, mpi_double_precision,
     &                my_rank - 1, 20, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine boundary_exchange_end(var)

      implicit none
      include '../par.for'
      include '../comm.par'
      include 'mpif.h'
      integer i, j
      integer status(MPI_Status_size)
      real*8 aux(inter*jmax), var(ptsx,jmax,msh)

      if (my_rank .lt. numproc) then
        ! Receiving the right columns from node + 1
        call MPI_Recv(aux, inter*jmax, mpi_double_precision,
     &                my_rank + 1, 249, MPI_COMM_WORLD, status, ierr)
        do i = 1, inter
          do j = 1, jmax
            var(i + ptsx - inter,j,1) = aux(j+(i-1)*jmax)
          end do
        end do
      end if
      if (my_rank .gt. 0) then
        ! Sending the left columns to node - 1
        do i = 1, inter
          do j = 1, jmax
            aux(j+(i-1)*jmax) = var(i + inter - inter + 1,j,1)
          end do
        end do
        call MPI_Send(aux, inter*jmax, mpi_double_precision,
     &                my_rank - 1, 249, MPI_COMM_WORLD, ierr)
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                          Poisson solver subroutine                           c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
