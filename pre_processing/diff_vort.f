ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c                 derivative calculations               c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program verifica

      implicit none
      include '../par.for'
      include '../comm.coef'
      integer i, j 
      real*8 stf_v, beta_fs_v, x, y
      real*8 dvx(imax,jmax), dvy(imax,jmax),
     &        dudy(imax,jmax), dvdx(imax,jmax),
     &        ux(imax,jmax), uy(imax,jmax),
     &        wz(imax,jmax)

      ! reads the boundary layer profile
      open(1,file='blasius.bin',form='unformatted')
      read(1) stf_v, beta_fs_v
      read(1) ux, uy, wz
      close(unit=1)

      ! verify if the files are compatible
      if(stf_v.ne.stf .or. beta_fs_v.ne.beta_fs) then
       write(*,*)
       write(*,*) 
     &  'Error! Binary files are incompatible with stf/FS parameters.'
       write(*,*)
       stop
      endif

      ! reads the derivative and Poisson coefficients
      open(1,file='coefs.bin',form='unformatted')
      read(1) stf_v
      read(1) fp_fd_coef
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
      read(1)
      read(1)
      read(1)
      read(1)
      close(unit=1)

      ! verify if the files are compatible
      if(stf_v.ne.stf) then
       write(*,*)
       write(*,*) 
     &  'Error! Binary files are incompatible with stf parameter.'
       write(*,*)
       stop
      endif

      ! mounts the lhs for the derivative calculation
      call derivs_k

      call dery(dudy,ux)
c     call derx(dvdx,uy)

      open (2, file = 'difference.dat',status = 'unknown')
      write(2,*) 'VARIABLES="x","y","vort1","vort2","diff"'
      write(2,*) 'ZONE I=',imax,', J=',jmax,', F=POINT'
 
      do j = 1, jmax
        if(stf.eq.1.d0) then
         y = dble(j-1) * dy
        else
         y = dy * (stf**(j-1)-1.d0)/(stf-1.d0)  
        endif        
        do i = 1, imax
          x = x0 + dble(i-1) * dx
          write(2,*)x,y,wz(i,j),dudy(i,j), 
     &             dabs(wz(i,j) - dudy(i,j))
        end do
      end do
      close(unit=2)


      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine derivs_k

      implicit none
      include '../par.for'
      real*8  a1x(imax),  b1x(imax),  c1x(imax)
      real*8  a1y(jmax),  b1y(jmax),  c1y(jmax)
      common/der1x/ a1x,b1x,c1x
      common/der1y/ a1y,b1y,c1y

      call coefx(a1x,b1x,c1x)
      call coef(a1y,b1y,c1y)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine derx(ddx,fc)

      ! first derivatives calculation in x direction
      implicit none
      include '../par.for'
      real*8 a1x(imax), b1x(imax), c1x(imax)
      real*8 fc(imax,jmax), ddx(imax,jmax)
      common/der1x/ a1x,b1x,c1x

      call rhsx(fc,ddx)
      call tridseqx(a1x,b1x,c1x,ddx)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dery(ddy,fc)

      ! first derivative calculation in y direction
      implicit none
      include '../par.for'
      real*8 a1y(jmax), b1y(jmax), c1y(jmax)
      real*8 fc(imax,jmax), ddy(imax,jmax)
      common/der1y/ a1y,b1y,c1y

      call rhsy(fc,ddy)
      call tridy(a1y,b1y,c1y,ddy)
      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsx(fc,rhs)
      
      ! RHS for the first derivative calculation in x direction
      implicit none
      include '../par.for'
      integer i, j, k
      real*8 rhs(imax,jmax), fc(imax,jmax)

      do j = 1, jmax
        rhs(1,j)=( - 74.d0 * fc(1,j) + 16.d0 * fc(2,j) +
     &                 72.d0 * fc(3,j) - 16.d0 * fc(4,j) +
     &                  2.d0 * fc(5,j) ) / ( 24.d0 * dx )
        
        rhs(2,j)=( - 406.d0 * fc(1,j) - 300.d0 * fc(2,j) +
     &                 760.d0 * fc(3,j) -  80.d0 * fc(4,j) +
     &                  30.d0 * fc(5,j) -   4.d0 * fc(6,j) ) /
     &             ( 120.d0 * dx )
        
        do i = 3, imax - 2
          rhs(i,j)=(           fc(i+2,j) - fc(i-2,j) +
     &                 28.d0 * ( fc(i+1,j) - fc(i-1,j) ) ) /
     &               ( 12.d0 * dx )
        end do
        
        rhs(imax-1,j)=( - 406.d0 * fc(imax,j)
     &                    - 300.d0 * fc(imax-1,j)
     &                    + 760.d0 * fc(imax-2,j)
     &                    -  80.d0 * fc(imax-3,j)
     &                    +  30.d0 * fc(imax-4,j)
     &                    -   4.d0 * fc(imax-5,j) ) /
     &                  ( - 120.d0 * dx )
        
        rhs(imax,j)=( - 74.d0 * fc(imax,j)
     &                  + 16.d0 * fc(imax-1,j)
     &                  + 72.d0 * fc(imax-2,j)
     &                  - 16.d0 * fc(imax-3,j)
     &                  +  2.d0 * fc(imax-4,j) ) /
     &                ( - 24.d0 * dx )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rhsy(fc,rhs)

      ! RHS for the first derivative calculation in y direction
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer i, j, k
      real*8 rhs(imax,jmax), fc(imax,jmax)
      real*8 ay

      do i = 1, imax
        ay = dy
        rhs(i,1) = (1.d0/ay) * ( fp_fd_coef(3) * fc(i,1) +
     &                             fp_fd_coef(4) * fc(i,2) +
     &                             fp_fd_coef(5) * fc(i,3) +
     &                             fp_fd_coef(6) * fc(i,4) +
     &                             fp_fd_coef(7) * fc(i,5) )
        
        rhs(i,2) = (1.d0/ay) * ( sp_fd_coef(4) * fc(i,1) + 
     &                             sp_fd_coef(5) * fc(i,2) +
     &                             sp_fd_coef(6) * fc(i,3) +  
     &                             sp_fd_coef(7) * fc(i,4) +
     &                             sp_fd_coef(8) * fc(i,5) +   
     &                             sp_fd_coef(9) * fc(i,6) )
        
       do j = 3, jmax - 2
         ay = dy * stf**(j-3)
         rhs(i,j) = (1.d0/ay) * ( cp_fd_coef(4) * fc(i,j-2) +
     &                              cp_fd_coef(5) * fc(i,j-1) +
     &                              cp_fd_coef(6) * fc(i,j)   +
     &                              cp_fd_coef(7) * fc(i,j+1) +
     &                              cp_fd_coef(8) * fc(i,j+2) )
       end do
       
       ay = dy * stf**(jmax-6)
       rhs(i,jmax-1) = (1.d0/ay) * ( pp_fd_coef(4)*fc(i,jmax)   +
     &                                 pp_fd_coef(5)*fc(i,jmax-1) +
     &                                 pp_fd_coef(6)*fc(i,jmax-2) +
     &                                 pp_fd_coef(7)*fc(i,jmax-3) +
     &                                 pp_fd_coef(8)*fc(i,jmax-4) +
     &                                 pp_fd_coef(9)*fc(i,jmax-5) )
       
       ay = dy * stf**(jmax-5)
       rhs(i,jmax) = (1.d0/ay) * ( lp_fd_coef(3) * fc(i,jmax)   +
     &                               lp_fd_coef(4) * fc(i,jmax-1) +
     &                               lp_fd_coef(5) * fc(i,jmax-2) +
     &                               lp_fd_coef(6) * fc(i,jmax-3) +
     &                               lp_fd_coef(7) * fc(i,jmax-4) )
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridy(a,b,c,rhs)

      ! solves tridiagonal matrix for the derivatives in y direction
      implicit none
      include '../par.for'
      integer i, j, k
      real*8 a(jmax), b(jmax), c(jmax), gam(jmax), bet
      real*8 rhs(imax,jmax), u(jmax)

      do i = 1, imax
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tridseqx(a,b,c,rhs)

      ! solves tridiagonal matrix for the derivatives in x direction
      implicit none
      include '../par.for'
      integer i, j, k
      real*8 a(imax), b(imax), c(imax), gam(imax), bet
      real*8 rhs(imax,jmax), u(imax)

      do j = 1, jmax
        bet  = b(1)
        u(1) = rhs(1,j) / bet
        do i = 2, imax
          gam(i) = c(i-1) / bet
          bet    = b(i) - a(i) * gam(i)
          u(i)   = ( rhs(i,j) - a(i) * u(i-1) ) / bet
        end do
        do i = imax - 1, 1, -1
          u(i) = u(i) - gam(i+1) * u(i+1)
        end do
        do i = 1, imax
          rhs(i,j) = u(i)
        end do
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coef(a,b,c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include '../par.for'
      include '../comm.coef'
      integer j
      real*8 a(jmax), b(jmax), c(jmax)

      a(1)      = 0.d0
      b(1)      = fp_fd_coef(1)
      c(1)      = fp_fd_coef(2)

      a(2)      = sp_fd_coef(1)
      b(2)      = sp_fd_coef(2)
      c(2)      = sp_fd_coef(3)

      do j = 3, jmax - 2
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine coefx(a,b,c)

      ! mount the LHS of the matrix for the first derivative
      implicit none
      include '../par.for'
      integer i
      real*8 a(imax), b(imax), c(imax)

      a(1)      = 0.d0
      b(1)      = 1.d0
      c(1)      = 4.d0

      a(2)      = 1.d0
      b(2)      = 6.d0
      c(2)      = 2.d0

      do i = 3, imax - 2
        a(i)    = 1.d0
        b(i)    = 3.d0
        c(i)    = 1.d0
      end do

      a(imax-1) = 2.d0
      b(imax-1) = 6.d0
      c(imax-1) = 1.d0

      a(imax)   = 4.d0
      b(imax)   = 1.d0
      c(imax)   = 0.d0

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c            end of derivative calculations             c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
