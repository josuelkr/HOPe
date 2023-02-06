cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                              c
c                   subroutines transforms binary data in asc                  c
c                                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program analysis 29052009

      implicit none
      include 'par.for'
      character*30 nome
      integer i, j, k, t, my_rank, shift, var, n, pos
      real*8  x, y, z,
     &        uxbs(ptsx,jmax), uybs(ptsx,jmax), wzbs(ptsx,jmax),
     &        thbs(ptsx,jmax),
     &        uxb(imax,jmax), uyb(imax,jmax), wzb(imax,jmax),
     &        thb(imax,jmax),
     &        uxp(imax,jmax,kphys),    wxp(imax,jmax,kphys),
     &        uyp(imax,jmax,kphys),    wyp(imax,jmax,kphys),
     &        uzp(imax,jmax,kphys),    wzp(imax,jmax,kphys),
     &        thp(imax,jmax,kphys),
     &        duxbdy(imax,jmax),       duxbdz(imax,jmax),
     &        dthbdy(imax,jmax),       dthbdz(imax,jmax),
     &        dthmdy(imax),            duxmdy(imax),
     &        duxdyp(imax,jmax,kphys), dthdyp(imax,jmax,kphys)
      real*8  t_duxmdy(imax,16), average_t_duxmdy(imax)
      real*8  t_dthmdy(imax,16), average_t_dthmdy(imax)
      real*8  lam, turb, ratio 
      complex*16 uxt(imax,jmax,kfour),    wxt(imax,jmax,kfour),
     &           uyt(imax,jmax,kfour),    wyt(imax,jmax,kfour),
     &           uzt(imax,jmax,kfour),    wzt(imax,jmax,kfour),
     &           tht(imax,jmax,kfour),   
     &            ux(ptsx,jmax,kfour),     wx(ptsx,jmax,kfour),
     &            uy(ptsx,jmax,kfour),     wy(ptsx,jmax,kfour),
     &            uz(ptsx,jmax,kfour),     wz(ptsx,jmax,kfour),
     &            th(ptsx,jmax,kfour),  duxdy(imax,jmax,kfour),
     &         dthdy(imax,jmax,kfour)

      call derivs_kt

      ! Boundary layer data
      do my_rank = 0, np - 1
        write(nome,'(a,i0.2,a)')'baseflow2D/based_',my_rank,'.bin'
        open(1, file = nome, form = 'unformatted')
        read(1) uxbs, uybs, wzbs, thbs
        close(unit = 1)
        shift = my_rank * (ptsx - inter - 1)
        do j = 1, jmax
          do i = 1, ptsx
            uxb(i+shift,j) = uxbs(i,j)
            thb(i+shift,j) = thbs(i,j)
          end do
        end do
      end do

      ! Disturbances variables data
      do n = 1, 16
        do my_rank = 0, np - 1
          write(nome,'(a,i0.2,a,i0.2)')'pert_',my_rank,'_',n
          open(2, file = nome, form = 'unformatted')
          read(2) ux, th
          close (unit = 2)
          shift = my_rank * (ptsx - inter - 1)
          do k = 1, kfour
            do j = 1, jmax
              do i = 1, ptsx
                uxt(i+shift,j,k) = ux(i,j,k)
                tht(i+shift,j,k) = th(i,j,k)
              end do
            end do
          end do
        end do

        call f_to_p(uxp, uxt, imax)
        call f_to_p(thp, tht, imax)
        
        call deryt(duxdy, uxt)
        call deryb(duxbdy, uxb)
        call deryt(dthdy, tht)
        call deryb(dthbdy, thb)
        
        call f_to_p(duxdyp, duxdy, imax)
        call f_to_p(dthdyp, dthdy, imax)
        
        do i = 1, imax
          duxmdy(i) = 0.d0
          dthmdy(i) = 0.d0
        end do
        do i = 1, imax
          do k = 1, kphys
            duxmdy(i) = duxmdy(i) + duxdyp(i,1,k)
            dthmdy(i) = dthmdy(i) + dthdyp(i,1,k)
          end do
          duxmdy(i) = duxmdy(i) / dble(kphys) + duxbdy(i,1)
          dthmdy(i) = dthmdy(i) / dble(kphys) + dthbdy(i,1)
          t_duxmdy(i,n) = duxmdy(i)
          t_dthmdy(i,n) = dthmdy(i)
        end do

        pos = (10.68d0-x0)/dx+1
        if ( n .eq .1 .or. n .eq. 5 .or. n .eq. 9 .or. n .eq. 13 ) then
          write(*,*) pos    
          write(nome,'(a,i0.2,a)')'crosscut_',n,'.dat'
          open (3, file = nome, status = 'unknown')
          write(3,*) 'VARIABLES="z","y","ux","th"'
          write(3,*) 'ZONE I=',kphys+1,', J=',jmax,', F=POINT'
        
          do j = 1, jmax
            if (stf .eq. 1.d0) then
              ! without stretching
              y = dble(j-1) * dy
             else
              ! with stretching
              y = dy * (stf**(j-1) - 1.d0) / (stf - 1.d0)
            end if
            do k = 1, kphys
              z = dble(k-1)*dz
              write(3,3) z, y, uxb(pos,j)+uxp(pos,j,k),
     &                         thb(pos,j)+thp(pos,j,k)
            end do
            write(3,3) z+dz, y, uxb(pos,j)+uxp(pos,j,1), 
     &                         thb(pos,j)+thp(pos,j,1)
          end do
          close (unit = 3)
        endif 
     
      end do

      average_t_duxmdy = sum(t_duxmdy, dim = 2) / 16d0
      average_t_dthmdy = sum(t_dthmdy, dim = 2) / 16d0


      open (4, file = 'cf.dat',status = 'unknown')
      write(4,*) 'VARIABLES=,"x","cf","laminar","turbulent","ratio"'
      write(4,*) 'ZONE T="cfs", I=',imax

      do i = 1, imax
        x = x0 + dble(i-1) * dx
        lam = 0.664d0*(U_1*x*L_1/N_1)**(-1.d0/2.d0)  
        turb = 0.027d0*(U_1*x*L_1/N_1)**(-1.d0/7.d0) 
        ratio = (average_t_dthmdy(i)*2d0/Re-lam)/(turb-lam)
        write(4,8) x, average_t_duxmdy(i)*2d0/Re, lam, turb, ratio  
      end do
      close(4)

      open (4, file = 'st.dat',status = 'unknown')
      write(4,*) 'VARIABLES=,"x","st","laminar","turbulent","ratio"'
      write(4,*) 'ZONE T="st", I=',imax

      do i = 1, imax
        x = x0 + dble(i-1) * dx
        lam = 0.3320d0*Pr**(-2.d0/3.d0)*(U_1*x*L_1/N_1)**(-1.d0/2.d0) 
        turb = 0.0296d0*Pr**(-2.d0/3.d0)*(U_1*x*L_1/N_1)**(-1.d0/5.d0)
        ratio = (average_t_dthmdy(i)/Re/Pr-lam)/(turb-lam)
        write(4,8) x, average_t_dthmdy(i)/Re/Pr, lam, turb, ratio  
      end do
      close(4)

      open (4, file = 'kappa.dat',status = 'unknown')
      write(4,*) 'VARIABLES=,"x","kappa"'
      write(4,*) 'ZONE T="kappa", I=',imax

      do i = 1, imax
        x = x0 + dble(i-1) * dx
c       write(4,2) x, - dtanh (3.0 * (x - 9.18045d0))   
c       write(4,2) x, 0.5d0 * (1- dtanh(3.0d0*(x-9.18045d0)))  
        write(4,2) x, 1d0  
      end do
      close(4)


    2 format(1x, 1d14.6, 1d17.9)
    3 format(1x, 2d14.6, 2d17.9)
    8 format(1x, 1d14.6, 4d17.9)
    5 format(1x, 1d14.6, 3d17.9)
    6 format(1x, 2d14.6, 6d17.9)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine derivs_kt

      ! calculate the tri-diagonal LHS for the derivative calculation
      implicit none
      include 'par.for'
      real*8  a1y(jmax), b1y(jmax), c1y(jmax)
      common/der1y/ a1y, b1y, c1y

      call coeft(a1y, b1y, c1y, jmax)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryt(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      complex*16 rhs(jmax), u(jmax),
     &           fc(imax,jmax,kfour), ddy(imax,jmax,kfour)
      common/der1y/ a1y, b1y, c1y

      do k = 1, kfour
        do i = 1, imax
          call rhsyt(fc, rhs, i, k)
          call tridyt(rhs, a1y, b1y, c1y, u)
          do j = 1, jmax
            ddy(i,j,k) = u(j)
          end do
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine deryb(ddy, fc)

      ! first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 a1y(jmax), b1y(jmax), c1y(jmax), rhs(jmax), u(jmax),
     &       fc(imax,jmax), ddy(imax,jmax)
      common/der1y/ a1y, b1y, c1y

      do i = 1, imax
        call rhsyb(fc, rhs, i)
        call tridyb(rhs, a1y, b1y, c1y, u)
        do j = 1, jmax
          ddy(i,j) = u(j)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyt(fc, rhs, i, k)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      complex*16 rhs(jmax), fc(imax,jmax,kfour)

      rhs(1) = ( - 74.d0 * fc(i,1,k) + 16.d0 * fc(i,2,k)
     &           + 72.d0 * fc(i,3,k) - 16.d0 * fc(i,4,k)
     &           +  2.d0 * fc(i,5,k) ) / ( 24.d0 * dy )

      rhs(2) = ( - 406.d0 * fc(i,1,k) - 300.d0 * fc(i,2,k)
     &           + 760.d0 * fc(i,3,k) -  80.d0 * fc(i,4,k)
     &           +  30.d0 * fc(i,5,k) -   4.d0 * fc(i,6,k) ) /
     &           ( 120.d0 * dy )

      do j = 3, jmax - 2
        rhs(j) = (           fc(i,j+2,k) - fc(i,j-2,k)
     &           + 28.d0 * ( fc(i,j+1,k) - fc(i,j-1,k) ) ) /
     &           ( 12.d0 * dy )
      end do

      rhs(jmax-1) = ( - 406.d0 * fc(i,jmax,k)
     &                - 300.d0 * fc(i,jmax-1,k)
     &                + 760.d0 * fc(i,jmax-2,k)
     &                -  80.d0 * fc(i,jmax-3,k)
     &                +  30.d0 * fc(i,jmax-4,k)
     &                -   4.d0 * fc(i,jmax-5,k) ) /
     &              ( - 120.d0 * dy )

      rhs(jmax) = ( - 74.d0 * fc(i,jmax,k)   + 16.d0 * fc(i,jmax-1,k)
     &              + 72.d0 * fc(i,jmax-2,k) - 16.d0 * fc(i,jmax-3,k)
     &              +  2.d0 * fc(i,jmax-4,k) ) / ( - 24.d0 * dy )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhsyb(fc, rhs, i)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include 'par.for'
      integer i, j, k
      real*8 rhs(jmax), fc(imax,jmax)

      rhs(1) = ( - 74.d0 * fc(i,1) + 16.d0 * fc(i,2)
     &           + 72.d0 * fc(i,3) - 16.d0 * fc(i,4)
     &           +  2.d0 * fc(i,5) ) / ( 24.d0 * dy )

      rhs(2) = ( - 406.d0 * fc(i,1) - 300.d0 * fc(i,2)
     &           + 760.d0 * fc(i,3) -  80.d0 * fc(i,4)
     &           +  30.d0 * fc(i,5) -   4.d0 * fc(i,6) ) /
     &           ( 120.d0 * dy )

      do j = 3, jmax - 2
        rhs(j) = (           fc(i,j+2) - fc(i,j-2)
     &           + 28.d0 * ( fc(i,j+1) - fc(i,j-1) ) ) /
     &           ( 12.d0 * dy )
      end do

      rhs(jmax-1) = ( - 406.d0 * fc(i,jmax)   - 300.d0 * fc(i,jmax-1)
     &                + 760.d0 * fc(i,jmax-2) -  80.d0 * fc(i,jmax-3)
     &                +  30.d0 * fc(i,jmax-4) -   4.d0 * fc(i,jmax-5) )
     &            / ( - 120.d0 * dy )

      rhs(jmax) = ( - 74.d0 * fc(i,jmax)   + 16.d0 * fc(i,jmax-1)
     &              + 72.d0 * fc(i,jmax-2) - 16.d0 * fc(i,jmax-3)
     &              +  2.d0 * fc(i,jmax-4) ) / ( - 24.d0 * dy )

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyt(rhs, a, b, c, u)

      ! solves the tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      complex*16 rhs(jmax), u(jmax)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
        u(j)   = ( rhs(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tridyb(rhs, a, b, c, u)

      ! solves the tridiagonal matrix for the derivatives in y direction
      implicit none
      include 'par.for'
      integer j
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet,
     &       rhs(jmax), u(jmax)

      bet  = b(1)
      u(1) = rhs(1) / bet
      do j = 2, jmax
        gam(j) = c(j-1) / bet
        bet    = b(j) - a(j) * gam(j)
        u(j)   = ( rhs(j) - a(j) * u(j-1) ) / bet
      end do
      do j = jmax - 1, 1, -1
        u(j) = u(j) - gam(j+1) * u(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine coeft(a, b, c, lmax)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      integer l, lmax
      real*8 a(lmax), b(lmax), c(lmax)
      
      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do l = 3, lmax - 2
        a(l)    = 1.d0
        b(l)    = 3.d0
        c(l)    = 1.d0
      end do

      a(lmax-1) = 2.d0
      b(lmax-1) = 6.d0
      c(lmax-1) = 1.d0

      a(lmax)   = 4.d0
      b(lmax)   = 1.d0
      c(lmax)   = 0.d0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
