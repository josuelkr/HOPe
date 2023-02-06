      program ose

      implicit none
      include 'par.f'
      integer jmax, jedge, j, i, imax
      complex*16 zeig(nm,jm), alfa, alfa0, re, freq, alfap, 
     &           z1, z2, z3, z4, dalfa
      real*8 eta(jm), h(jm), hp(jm), y(jm), yp(jm), x, reout,
     &       wslope, etamax, bl, re0, w0, beta, w0in, betain, rein, 
     &       fpp(10000), dx, Re_dns, beta_dns, omega_dns, beta_fsin, 
     &       betafs(10000), delta_re, beta_fs, wslopein, rea
      common /data/ re, alfa, alfa0, freq, beta, re0, w0, beta_fs,
     &       /blas/ wslope, etamax,
     &       /cons/ z1, z2, z3, z4
      namelist/osedata/ re0, w0, alfa0, jedge, bl, etamax, beta, wslope

c************************************************************************
c     this code solves the Orr-Sommerfeld equation using v (normal
c     velocity) as dependent variable.
c     the code can use either rkdumb or odeint. For that
c     change the call statment on init.f, and change the 
c     assignment of yp=eta in the write statement bellow.
c************************************************************************

      open(unit = 1, file = 'ose3d.in', status = 'old')
        read(1, nml = osedata)
      close(unit = 1)

      open(2, file = 'datalst.dat', form = 'formatted')
        read(2,*) Re_dns, omega_dns, beta_dns, imax, dx
      close(unit = 2)

      open(3, file = 'fpp.dat', form = 'formatted')
      open(4, file = '../beta_fs.dist', form = 'formatted')
      do i = 1, imax
        read(3,*) fpp(i)
        read(4,*) betafs(i)
      end do
      close(unit = 3)
      close(unit = 4)

      w0in    = w0
      betain  = beta
      rein    = re0
      beta_fs = 0.d0

      call geom(eta, y, h, hp, bl, etamax, jmax, jedge)

      ! loop para Reynolds
      do i = 1, 1000
        re0  = rein + dble(i-1) * (re_dns - rein) / 1000.d0
        call evalue (zeig, jmax, eta, yp)
      end do

      ! loop para omega
      do i = 1, 1000
        w0   = w0in + dble(i-1) * (omega_dns - w0in) / 1000.d0
        call evalue (zeig, jmax, eta, yp)
      end do

      ! loop para beta (spanwise wavenumber)
      do i = 1, 1000
        beta = betain + dble(i-1) * (beta_dns - betain) / 1000.d0
        call evalue (zeig, jmax, eta, yp)
      end do
    
      ! loop para fpp e beta_fs
      wslopein  = wslope
      beta_fsin = 0.d0
      do i = 1, 1000
        wslope  = wslopein  + dble(i-1) * (fpp(1)    - wslopein)  / 1000.d0
        beta_fs = beta_fsin + dble(i-1) * (betafs(1) - beta_fsin) / 1000.d0
        call evalue (zeig, jmax, eta, yp)
      end do

      rein   = re_dns
      re0    = re_dns
      w0in   = w0
      betain = beta

      open (1, file = 'lst.dat', status = 'unknown')
      write(1,*) 'VARIABLES="x","alfa"'
      write(1,*) 'ZONE T=LST, I=',imax

      do i = 1, imax
        x       = 1.d0 + dble(i-1) * dx
        wslope  = fpp(i)
        beta_fs = betafs(i) 
        re0     = rein * dsqrt(x) 
        w0      = w0in * re0 / rein
        beta    = betain * re0 / rein
        call evalue (zeig, jmax, eta, yp)
c       call eigen (eta, y, h, hp, zeig, jmax, jedge)
        x       = (re0 / rein)**2
        write(*,*) i, x, re0, alfa0
        write(1,3) x, -dimag(alfa0)*rein**2/re0
      end do
      close (unit = 1)
    3 format(1x, 2e17.9)

      stop
      end

**********************************************************************      

      subroutine derivs(x, zy, zd)

      implicit none
      include 'par.f'
      complex*16 utot, uptot, upptot, zy(nm), zd(nm),
     &           re, alfa, alfa0, freq,
     &           z1, z2, z3, z4, z5, z6, z7, z8, z9
      real*8 beta, re0, w0, x, m, beta_fs
      common /data/ re, alfa, alfa0, freq, beta, re0, w0, beta_fs,
     &       /cons/ z1, z2, z3, z4

      m = beta_fs / (2.d0 - beta_fs)

      zd(1)  = zy(2)
      zd(2)  = zy(3)
      zd(3)  = - dcmplx((m + 1.d0)/ 2.d0) * zy(1) * zy(3)
     &         - dcmplx(m) * ( dcmplx(1.d0) - zy(2)**2 ) 
      utot   = zy(2) 
      uptot  = zy(3)
      upptot = zd(3) 

c     z1     = re * i * freq
c     z2     = alfa**2
c     z3     = re * i * alfa
c     z4     = z2 * (- z2 + z1)
      z5     = + z3 * utot
      z6     = - z1 + z2 * dcmplx(2.d0) + z5
      z7     = + z4 - z3 * (upptot + z2 * utot)
      z8     = - z1 + z2 + z5
      z9     = + dcmplx(re0 * beta) * uptot
      
      zd(4)  = zy(5)
      zd(5)  = zy(6)
      zd(6)  = zy(7)
      zd(7)  = z6 * zy(6) + z7 * zy(4)
      zd(8)  = zy(9)
      zd(9)  = z8 * zy(8) + z9 * zy(4) 

      zd(10) = zy(11)
      zd(11) = zy(12)
      zd(12) = zy(13)
      zd(13) = z6 * zy(12) + z7 * zy(10)
      zd(14) = zy(15)
      zd(15) = z8 * zy(14) + z9 * zy(10) 

      zd(16) = zy(17)
      zd(17) = zy(18)
      zd(18) = zy(19)
      zd(19) = z6 * zy(18) + z7 * zy(16)
      zd(20) = zy(21)
      zd(21) = z8 * zy(20) + z9 * zy(16) 

      return
      end


**********************************************************************      

      subroutine deter(zy, zdet)

      implicit none
      integer j, indx(6)
      include 'par.f'
      complex*16 alfa, alfa0, gamma, re, freq,
     &           za(6,6), zy(nm), zdet, z1, z2, z3, z4
      real*8 wslope, etamax, rea, aim, re0, beta, w0
      common /data/ re, alfa, alfa0, freq, beta, re0, w0,
     &       /blas/ wslope, etamax,
     &       /cons/ z1, z2, z3, z4

      gamma = (freq - alfa) * im * re
      gamma = zsqrt (z2 - gamma)
      rea   = dreal(gamma)
      aim   = dimag(gamma)

      if (rea .lt. 0.d0) then          !check if real parte of gamma is > 0
        rea = - rea 
      endif
      if ((rea .le. 0.d0) .and. (aim .gt. 0.d0)) then
        aim = -aim
      endif
      gamma = dcmplx(rea,aim)

      za(1,1) = zy(4)
      za(2,1) = zy(5)
      za(3,1) = zy(6)
      za(4,1) = zy(7)
      za(5,1) = zy(8)
      za(6,1) = zy(9)

      za(1,2) = zy(10)
      za(2,2) = zy(11)
      za(3,2) = zy(12)
      za(4,2) = zy(13)
      za(5,2) = zy(14)
      za(6,2) = zy(15)

      za(1,3) = zy(16)
      za(2,3) = zy(17)
      za(3,3) = zy(18)
      za(4,3) = zy(19)
      za(5,3) = zy(20)
      za(6,3) = zy(21)

      za(1,4) = - dcmplx(1.d0)
      za(2,4) =   alfa
      za(3,4) = - alfa * alfa
      za(4,4) =   alfa * alfa * alfa
      za(5,4) =   dcmplx(0.d0)
      za(6,4) =   dcmplx(0.d0)

      za(1,5) = - dcmplx(1.d0)
      za(2,5) =   gamma
      za(3,5) = - gamma * gamma
      za(4,5) =   gamma * gamma * gamma
      za(5,5) =   dcmplx(0.d0)
      za(6,5) =   dcmplx(0.d0)

      za(1,6) =   dcmplx(0.d0)
      za(2,6) =   dcmplx(0.d0)
      za(3,6) =   dcmplx(0.d0)
      za(4,6) =   dcmplx(0.d0)
      za(5,6) = - dcmplx(1.d0)
      za(6,6) =   gamma

      call ludcmp(za, 6, 6, indx, zdet)

      do j = 1, 6
        zdet = zdet * za(j,j)
      end do

      return
      end


**********************************************************************      

      subroutine eigen(eta, y, h, hp, zeig, jmax, jedge)

      implicit none
      integer jmax, jedge, icont, i, j, knorm, iflag, jcut
      include 'par.f'
      complex*16 zeig(nm,jm), za(6,6), zb(6), z1(6), re, freq,
     &           zbt(6), zd, zfi(6,jm), alfa, alfa0, gamma,
     &           expoa, expog, vmax
      real*8 eta(jm), y(jm), h(jm), hp(jm), fo(3,jm), re0, beta,
     &       wslope, etamax, anorm, bmax, bmag, c1, test, w0
      integer indx(6)
      common /data/ re, alfa, alfa0, freq, beta, re0, w0,
     &       /blas/ wslope, etamax

c************************************************************************
c    inverse iteration subroutine.
c    In this subroutine the outer edge boundary condition is
c    matched with the three orthonormal solutions in order to
c    calculate the eigenfunction for a given Re, freq and alfa.
c    The system of equations is given by:
c
c    A*fI + B*fII + C*fIII - D*exp(-alfa*etamax) -
c         - E*exp(-gamma*etamax) -F*dexp(-gamma*etamax) = 0.
c
c    and respective derivatives.
c
c    fI, fII and fIII are the three solution vectors.
c    the assymptotic solutions are:
c
c            v = D * exp(-alfa*etamax) + E * exp(-gamma*etamax)
c            w = F * exp(-gamma*etamax)
c
c************************************************************************

      call invit(zb, zeig, jmax, gamma)

c..............................................fo is the blasius function
c.............................zfi(1) to zfi(4) are v, and its derivatives
c..................zfi(5) and zfi(6) are the vorticity and its derivative

c......note..............................................................
c.................out of the rk4/ortho routines the profiles have a small
c...............................discontinuity. To ajust for that zeig(20)
c...............................and zeig(21) are disconsidered since they
c.......................................are supposed to go to zero anyway
c....................................also zb(5) is set to zero to ajust a 
c............................................discontinuity on Im(zfi(4)).
c.....................................this discontinuity is there because
c........................................zb's are being calculated with a
c.............................discontinuous profile at the outer boundary
c..............................search for "!" to locate the modifications
c........................................................................

      iflag = jmax
      do j = 1, jmax
        zfi(1,j) = zb(1) * zeig(4,j) + zb(2) * zeig(10,j)
     &           + zb(3) * zeig(16,j)
        zfi(2,j) = zb(1) * zeig(5,j) + zb(2) * zeig(11,j)
     &           + zb(3) * zeig(17,j)
        zfi(3,j) = zb(1) * zeig(6,j) + zb(2) * zeig(12,j)
     &           + zb(3) * zeig(18,j)
        zfi(4,j) = zb(1) * zeig(7,j) + zb(2) * zeig(13,j)
     &           + zb(3) * zeig(19,j)
        zfi(5,j) = zb(1) * zeig(8,j) + zb(2) * zeig(14,j)
!    &           + zb(3) * zeig(20,j)
        zfi(6,j) = zb(1) * zeig(9,j) + zb(2) * zeig(15,j)
!    &           + zb(3) * zeig(21,j)
        fo(1,j)  = dreal(zeig(1,j))
        fo(2,j)  = dreal(zeig(2,j))
        fo(3,j)  = dreal(zeig(3,j))
        if(fo(2,j) .gt. 1.d0) fo(2,j) = 1.d0
        if((fo(3,j) .lt. 1.d-20) .and. (j .gt. jmax/2)) fo(3,j) = 0.d0
        if(zabs(zfi(1,j)) .gt. zabs(vmax)) vmax = zfi(1,j)
      end do

c..................................extrapolate the Blasius Boundary Layer
      c1 = eta(jmax) - fo(1,jmax)
      do j = jmax + 1, jedge - 1
        fo(1,j) = eta(j) - c1
        fo(2,j) = 1.d0
        fo(3,j) = 0.d0
      end do

c....................................extrapolate the eigenfunctions using
c.......................................................exponential decay
      do j = jmax + 1, jedge - 1
        expoa    = zexp(- alfa  * dcmplx(eta(j)) )
        expog    = zexp(- gamma * dcmplx(eta(j)) )
        zfi(1,j) = + zb(4) * expoa
     &             + zb(5) * expog
        zfi(2,j) = - zb(4) * expoa * alfa 
     &             - zb(5) * expog * gamma 
        zfi(3,j) = + zb(4) * expoa * alfa * alfa
     &             + zb(5) * expog * gamma * gamma
        zfi(4,j) = - zb(4) * expoa * alfa * alfa * alfa
!    &             - zb(5) * expog * gamma * gamma * gamma   ! ignore
        zfi(5,j) = dcmplx(0.d0) ! + zb(6) * expog            ! ignore
        zfi(6,j) = dcmplx(0.d0) ! - zb(6) * expog * gamma    ! ignore
      end do

      do j = 1, jedge - 1
        do i = 1, 6
          zfi(i,j) = zfi(i,j) / vmax
        end do
      end do

c..............................................write the Blasius solution
c     open(unit = 3, file = 'mean.dat', status = 'unknown')
c     do j = 1, jedge - 1
c        write(3,*) fo(1,j), fo(2,j), fo(3,j)
c         write(99,*) real(eta(j)), real(fo(1,j)), real(fo(2,j)), 
c    &    real(fo(3,j))
c     end do
c     close(unit = 3)

      call out(eta, y, zfi, fo, h, hp, jedge)

      return      
      end

**********************************************************************      

      subroutine evalue(zeig, jmax, eta, yp)

      implicit none
      integer icount, jmax, j
      include 'par.f'
      complex*16 zeig(nm,jm), alfa, alfa0, alfa2, re, freq,
     &           zx1, zfl, zx2, zxl, zf, zrtsec,
     &           zswap, zdx, zdet1, zdet2,
     &           z1, z2, z3, z4
      real*8 eta(jm), yp(jm), ratr, rati, beta, re0, w0
      common /data/ re, alfa, alfa0, freq, beta, re0, w0,
     &       /cons/ z1, z2, z3, z4

      icount = 0
      alfa   = zsqrt( alfa0**2 + dcmplx(beta**2) )
      alfa2  = dcmplx(1.0001d0) * alfa
      re     = re0 * alfa0 / alfa
      freq   = w0  * alfa  / alfa0
      z1     = re * im * freq
      z2     = alfa * alfa
      z3     = im * re * alfa
      z4     = z2 * (- z2 + z1)

      call init(zdet1, zeig, jmax, eta, yp)
c     write(*,*) alfa0

      zx1   = alfa
      zfl   = zdet1
      alfa  = alfa2
      alfa0 = zsqrt(alfa**2 - dcmplx(beta**2))
      re    = re0 * alfa0 / alfa
      freq  = w0  * alfa  / alfa0
      z1    = re * im * freq
      z2    = alfa * alfa
      z3    = im * re * alfa
      z4    = z2 * (- z2 + z1)

      call init(zdet2, zeig, jmax, eta, yp)
c     write(*,*) alfa0

      zx2 = alfa
      zf  = zdet2
      if (zabs(zfl) .lt. zabs(zf)) then
        zrtsec = zx1
        zxl    = zx2
        zswap  = zfl
        zfl    = zf
        zf     = zswap
       else
        zxl    = zx1
        zrtsec = zx2
      end if

      do j = 1, 100
        zdx    = (zxl - zrtsec) * zf / (zf - zfl)
        zxl    = zrtsec
        zfl    = zf
        zrtsec = zrtsec + zdx

        alfa   = zrtsec
        alfa0  = zsqrt(alfa**2 - dcmplx(beta**2))
        re     = re0 * alfa0 / alfa
        freq   = w0  * alfa  / alfa0
        z1     = re * im * freq
        z2     = alfa * alfa
        z3     = alfa * im * re
        z4     = z2 * (- z2 + z1)

        call init(zf, zeig, jmax, eta, yp)
c       write(*,*) alfa0

        ratr = dreal(zdx) / dreal(zrtsec)
        if (dabs(dimag(zrtsec)) .ge. 1.d-6) then
          rati = dimag(zdx) / dimag(zrtsec)
          if ((dabs(ratr) .le. 5.d-5) .and.
     &        (dabs(rati) .le. 1.d-2)) return
         else
          rati = dimag(zdx)
          if ((dabs(ratr) .le. 5.d-5) .and.
     &        (dabs(rati) .le. 5.d-5)) return
        end if
      end do

      write(*,*) 'too many iterations in evalue'
      stop

      return
      end

**********************************************************************      

      subroutine geom(eta, y, h, hp, bl, etamax, jmax, jedge)

      implicit none
      integer jmax, jedge, j
      include 'par.f'
      real*8 y(jm), eta(jm), h(jm), hp(jm), dy, etal, bl, etamax

c...............................................for mapping 0-infinty use
      jmax = idint(dble(jedge) * bl)
      dy   = 1.d0 / dble(jedge-1)
      etal = etamax * (1.d0 - dble(jmax-1) * dy) / dble(jmax-1) / dy

      do j = 1, jedge - 1
        y(j)   =   dble(j-1) * dy
        eta(j) =   etal * y(j) / (1.d0 - y(j))
        h(j)   =   (eta(j) + etal)**2 / etal
        hp(j)  = - 2.d0 * etal / (eta(j) + etal)**3
      end do
      y(jedge)  = 1.d0

c     write(*,*)
c     write(*,*) 'L      = ',etal
c     write(*,*) 'jedge  = ',jedge
c     write(*,*) 'jmax   = ',jmax
c     write(*,*) 'etamax = ',etamax, eta(jmax)
c     write(*,*)

      eta(jmax) = etamax

c     open(unit = 21, file = 'eta.dat', status = 'unknown')
c     open(unit = 22, file = 'h.dat', status = 'unknown')
c     do j = 1, jedge-1
c        write(21,*) y(j), eta(j)
c        write(22,*) h(j), hp(j)
c     end do
c     close(unit = 21)
c     close(unit = 22)

      return
      end

**********************************************************************      

      subroutine init (zdet, zeig, jmax, eta, yp)

      implicit none
      integer nvar, jmax
      include 'par.f'
      complex*16 zy(nm), zystart(nm), zeig(nm,jm), zdet
      real*8 eta(jm), yp(jm), x1, x2, eps, h1, hmin

      nvar = nm
      x1   = eta(1)
      x2   = eta(jmax)
      eps  = 1.d-4
      h1   = (eta(2) - eta(1))
      hmin = (eta(2) - eta(1)) * 0.25d0

      call start  (zystart)
      call rkdumb (zystart, nvar, zy, zeig, jmax, eta)
c     call odeint (zystart, nvar, x1, x2, eps, h1, hmin, zy, zeig, jmax, yp)

      call deter (zy, zdet)

      return
      end

**********************************************************************      

      subroutine invit(zb, zeig, jmax, gamma)

      implicit none
      integer jmax, icont, j, knorm, i, ii
      include 'par.f'
      complex*16 zeig(nm,jm), za(6,6), zb(6), z1(6), re, freq,
     &           zbt(6), zd, alfa, alfa0, gamma, z2
      real*8 re0, beta, w0, wslope, etamax, 
     &       anorm, bmax, bmag, test
      integer indx(6)
      common /data/ re, alfa, alfa0, freq, beta, re0, w0,
     &       /blas/ wslope, etamax

c************************************************************************
c    inverse iteration subroutine.
c    In this subroutine we calculate coeficients defined by the matching
c    of the solution at the outer edge to a exponential decay solution
c    
c    The system of equations is given by:
c
c    A*fI + B*fII + C*fIII + D*exp(-alfa*etamax) +
c    E*exp(-gamma*etamax) + F*dexp(-gamma*etamax) = 0.
c
c    and respective derivatives.
c    The subroutine uses inverse iteration to find the 
c    constants A, B, C, D, E, and F
c
c    fI, fII and fIII are the three solution vectors.
c    the assymptotic solutions are:
c
c    v = D * exp(-alfa * etamax) + E * exp(-gamma * etamax)
c    w = F * exp(-gamma * etamax)
c
c************************************************************************

      gamma   = alfa * alfa + (alfa - freq) * im * re
      gamma   = zsqrt(gamma)

      za(1,1) =   zeig(4,jmax)
      za(2,1) =   zeig(5,jmax) 
      za(3,1) =   zeig(6,jmax) 
      za(4,1) =   zeig(7,jmax)
      za(5,1) =   zeig(8,jmax)
      za(6,1) =   zeig(9,jmax)

      za(1,2) =   zeig(10,jmax)
      za(2,2) =   zeig(11,jmax)
      za(3,2) =   zeig(12,jmax) 
      za(4,2) =   zeig(13,jmax)
      za(5,2) =   zeig(14,jmax)
      za(6,2) =   zeig(15,jmax)

      za(1,3) =   zeig(16,jmax)
      za(2,3) =   zeig(17,jmax) 
      za(3,3) =   zeig(18,jmax) 
      za(4,3) =   zeig(19,jmax)
      za(5,3) =   zeig(20,jmax)
      za(6,3) =   zeig(21,jmax)

      za(1,4) = - dcmplx(1.d0)       * zexp(- alfa * dcmplx(etamax))
      za(2,4) =   alfa               * zexp(- alfa * dcmplx(etamax))
      za(3,4) = - alfa * alfa        * zexp(- alfa * dcmplx(etamax))
      za(4,4) =   alfa * alfa * alfa * zexp(- alfa * dcmplx(etamax))
      za(5,4) =   dcmplx(0.d0)
      za(6,4) =   dcmplx(0.d0)

      za(1,5) = - dcmplx(1.d0)          * zexp(- gamma * dcmplx(etamax))
      za(2,5) =   gamma                 * zexp(- gamma * dcmplx(etamax))
      za(3,5) = - gamma * gamma         * zexp(- gamma * dcmplx(etamax))
      za(4,5) =   gamma * gamma * gamma * zexp(- gamma * dcmplx(etamax))
      za(5,5) =   dcmplx(0.d0)
      za(6,5) =   dcmplx(0.d0)

      za(1,6) =   dcmplx(0.d0)
      za(2,6) =   dcmplx(0.d0)
      za(3,6) =   dcmplx(0.d0)
      za(4,6) =   dcmplx(0.d0)
      za(5,6) = - dcmplx(1.d0) * zexp(- gamma * dcmplx(etamax))
      za(6,6) =   gamma        * zexp(- gamma * dcmplx(etamax))

c.......................the method of solution for the homogenios systems 
c........................is inverse iteration where a initial guess for b 
c......................in Ax=b is updted with the solution x, normalised.

      zb(1)   = (1.d0,0.d0) 
      zb(2)   = (1.d-3,0.d0)
      zb(3)   = (1.d-6,0.d0)
      zb(4)   = (1.d0,0.d0)
      zb(5)   = (1.d-6,0.d0)
      zb(6)   = (1.d-6,0.d0)

      call ludcmp(za, 6, 6, indx, zd)

      icont = 0

      test = 1.d0
      do while (test .ge. 1.d-14)
        do j = 1, 6
          zbt(j) = zb(j)
        end do
        icont = icont + 1

        call lubksb(za, 6, 6, indx, zb)

        anorm = 0.d0
        do j = 1, 6
          anorm = anorm + zb(j) * dconjg(zb(j))
        end do

        anorm = dsqrt(anorm)
        do j = 1, 6
          zb(j) = zb(j) / dcmplx(anorm)
        end do

c..................................test for convergency of the solution. 
c.......(the norm of the diff between the present and previous solution)

        bmax = 0.d0
        do j = 1, 6
          bmag = zb(j) * dconjg(zb(j))
          bmag = dsqrt(bmag)
          if (bmag .gt. bmax)then
            knorm = j
            bmax  = bmag
          end if
        end do

        test = 0.d0
        do j = 1, 6
          z1(j) = zb(j) / zb(knorm) - zbt(j) / zbt(knorm)
          test = test + z1(j) * dconjg(z1(j))
        end do
        test = dsqrt(test)

        if(icont .ge. 100)then
          write(*,*)
          write(*,*) 'to many iterations to find the eigenfunctions'
          write(*,*) icont, 'iterations'
          write(*,*) test, 'maximum accuracy level'
          return
        endif
      end do           ! do while (test .ge. 1.d-14)

c     write(*,*)
c     write(*,*) 'convergency is accurate to ',test
c     write(*,*) 'after ',icont,' iterations'

      return
      end

**********************************************************************      

      subroutine lubksb(a, n, np, indx, b)

      implicit none
      integer n, np, ii, ll, j, i, indx(n)
      complex*16 a(np,np), b(np), sum

      ii = 0
      do i = 1, n
        ll    = indx(i)
        sum   = b(ll)
        b(ll) = b(i)
        if (ii .ne. 0) then
          do j = ii, i - 1
            sum = sum - a(i,j) * b(j)
          end do
         else if (sum .ne. dcmplx(0.d0)) then
          ii = i
        endif 
        b(i) = sum
      end do

      do i = n, 1, -1
        sum = b(i)
        do j = i + 1, n
          sum = sum - a(i,j) * b(j)
        end do
        b(i) = sum / a(i,i)
      end do

      return
      end

**********************************************************************      

      subroutine ludcmp (a, n, np, indx, d)

      implicit none
      integer n, np, indx(n), i, j, k, imax
      complex*16 a(np,np), d, sum, dummy
      real*8 vv(21), aamax, dum, tiny
      parameter (tiny = 1.d-10)

      d = dcmplx(1.d0)
      do i = 1, n
        aamax = 0.d0
        do j = 1, n
          if (zabs(a(i,j)) .gt. aamax) aamax = zabs(a(i,j))
        end do
        if (aamax .eq. 0.d0) then
          write(*,*) 'singular matrix'
          stop
        end if
        vv(i) = 1.d0 / aamax
      end do

      do j = 1, n
        do i = 1, j - 1
          sum = a(i,j)
          do k = 1, i - 1
            sum = sum - a(i,k) * a(k,j)
          end do
          a(i,j) = sum
        end do

        aamax = 0.d0
        do i = j, n
          sum = a(i,j)
          do k = 1, j - 1
            sum = sum - a(i,k) * a(k,j)
          end do
          a(i,j) = sum
          dum    = vv(i) * zabs(sum)
          if (dum .ge. aamax) then
            imax  = i
            aamax = dum
          end if
        end do

        if (j .ne. imax) then
          do k = 1, n
            dummy     = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k)    = dummy
          end do
          d = - d
          vv(imax) = vv(j)
        end if

        indx(j) = imax
        if (a(j,j) .eq. dcmplx(0.d0)) a(j,j) = tiny
        if (j .ne. n) then
          dummy = dcmplx(1.d0) / a(j,j)
          do i = j + 1, n
            a(i,j) = a(i,j) * dummy
          end do
        end if
      end do

      return
      end

**********************************************************************      

      subroutine odeint(zystart, nvar, x1, x2, eps, h1, hmin, zy,
     &                  zeig, jmax, yp)

      implicit none
      integer nvar, jmax, maxstp, i, j, jcount, jj
      include 'par.f'
      complex*16 zystart(nvar), zyscal(nm), zy(nm), 
     &           zdydx(nm), zeig(nm,jm)
      real*8 x, x1, x2, eps, h, h1, hmin, hnext, xsav, tiny, yp(jm)
      parameter (tiny = 1.d-30, maxstp = jm)

      x      = x1
      h      = dsign(h1,x2-x1)
      jcount = 1
      yp(1)  = x
 
      do i = 1, nvar
        zy(i) = zystart(i)
      end do

      do j = 1, maxstp
        call derivs(x, zy, zdydx)
        do i = 1, nvar
          zyscal(i) = zabs(zy(i))
     &              + zabs(dcmplx(h) * zdydx(i) + dcmplx(tiny))
        end do

        if ((x+h-x2)*(x+h-x1) .gt. 0.d0) h = x2 - x

        call rkqs(zy, zdydx, nvar, x, h, eps, zyscal, hnext)

        call ortho(zy, zeig, j+1, jmax)

        yp(j+1) = x
        jcount  = j + 1
        do i = 1, nvar
          zeig(i,j+1) = zy(i)
        end do

        if ((x-x2)*(x2-x1) .ge. 0.d0) then   ! check if it is done
          do i = 1, nvar
            zystart(i) = zy(i)
          end do
c         write(*,*) jcount
          return
        end if

        if (abs(hnext) .lt. hmin) hnext = hmin

        h = hnext
      end do

      write(*,*) 'Too many steps in odeint'
      stop
      end 

**********************************************************************      

      subroutine rkqs(y, dydx, n, x, htry, eps, yscal, hnext)

      implicit none
      integer n, nmax, i
      parameter (nmax = 21)
      real*8 eps, hnext, htry, x,
     &       errmax, h, xnew, safety, pgrow, pshrnk, errcon
      complex*16 dydx(n), y(n), yscal(n), yerr(nmax), ytemp(nmax) 
      parameter (safety = 0.9d0, pgrow = - 0.2d0, pshrnk = - 0.25d0,
     &           errcon = 1.89d-4)

      h = htry

 10   call rkck(y, dydx, n, x, h, ytemp, yerr)

      errmax = 0.d0
      do i = 1, n
        errmax = max(errmax, zabs(yerr(i) / yscal(i)))
      end do

      errmax = errmax / eps
      if(errmax .gt. 1.) then
        h = safety * h * (errmax**pshrnk)
        if (h .lt. 0.1*h) then
          h = .1d0 * h
        end if
        xnew = x + h
        if (xnew .eq. x) pause 'stepsize underflow in rkqs'
        go to 10
       else
        if (errmax .gt. errcon) then
          hnext = safety * h * (errmax**pgrow)
         else
          hnext = 5.d0 * h
        end if
        x = x + h
        do i = 1, n
          y(i) = ytemp(i)
        end do
        return
      end if

      return
      end

**********************************************************************      

      subroutine rkck(y, dydx, n, x, h, yout, yerr)

      implicit none
      integer n, nmax, i
      parameter (nmax = 21)
      real*8 h, x,
     &       a2, a3, a4, a5, a6, b21, b31, b32, b41, b42, b43, b51,
     &       b52, b53, b54, b61, b62, b63, b64, b65, c1, c3, c4, 
     &       c6, dc1, dc3, dc4, dc5, dc6
      complex*16 dydx(n), y(n), yerr(n), yout(n), ytemp(nmax),
     &           ak2(nmax), ak3(nmax), ak4(nmax), 
     &           ak5(nmax), ak6(nmax)
      parameter (a2 = 0.2d0, a3 = 0.3d0, a4 = 0.6d0, a5 = 1.d0,
     &           a6 = 0.875d0, b21 = 0.2d0, b31 = 3.d0 / 40.d0,
     &           b32 = 9.d0 / 40.d0, b41 = 0.3d0, b42 = - 0.9d0,
     &           b43 = 1.2d0, b51 = - 11.d0 / 54.d0, b52 = 2.5d0,
     &           b53 = - 70.d0 / 27.d0, b54 = 35.d0 / 27.d0,
     &           b61 = 1631.d0 / 55296.d0, b62 = 175.d0 / 512.d0, 
     &           b63 = 575.d0 / 13824.d0, b64 = 44275.d0 / 110592.d0,
     &           b65 = 253.d0 / 4096.d0, c1 = 37.d0 / 378.d0,
     &           c3 = 250.d0 / 621.d0, c4 = 125.d0 / 594.d0,
     &           c6 = 512.d0 / 1771.d0, dc1 = c1 - 2825.d0 / 27648.d0,
     &           dc3 = c3 - 18575.d0 / 48384.d0,
     &           dc4 = c4 - 13525.d0 / 55296.d0,
     &           dc5 = - 277.d0 / 14366.d0, dc6 = c6 - .25d0)

      do i = 1, n
        ytemp(i) = y(i) + b21 * h * dydx(i)
      end do
      
      call derivs(x+a2*h, ytemp, ak2)
      
      do i = 1, n
        ytemp(i) = y(i) + h * (b31 * dydx(i) + b32 * ak2(i))
      end do

      call derivs(x+a3*h,ytemp,ak3)

      do i=1,n
        ytemp(i) = y(i) + h * ( b41 * dydx(i) + b42 * ak2(i)
     &                        + b43 * ak3(i) )
      end do

      call derivs(x+a4*h, ytemp, ak4)

      do i = 1, n
        ytemp(i) = y(i) + h * ( b51 * dydx(i) + b52 * ak2(i)
     &                        + b53 * ak3(i)  + b54 * ak4(i) )
      end do

      call derivs(x+a5*h, ytemp, ak5)

      do i = 1, n
        ytemp(i) = y(i) + h * ( b61 * dydx(i) + b62 * ak2(i)
     &                        + b63 * ak3(i)  + b64 * ak4(i)
     &                        + b65 * ak5(i) )
      end do

      call derivs(x+a6*h, ytemp, ak6)

      do i = 1, n
        yout(i) = y(i) + h * ( c1 * dydx(i) + c3 * ak3(i)
     &                       + c4 * ak4(i)  + c6 * ak6(i) )
      end do

      do i = 1, n
        yerr(i) = h * ( dc1 * dydx(i) + dc3 * ak3(i) + dc4 * ak4(i)
     &                + dc5 * ak5(i)  + dc6 * ak6(i) )
      end do

      return
      end

**********************************************************************      

      subroutine ortho(zy, zeig, jj, jmax)

      implicit none
      integer jj, i, j, k, jmax
      include 'par.f'
      complex*16 zy(nm), zeig(nm,jm), gamma12, gamma13, gamma23
      real*8 amag1, amag2, amag3

c************************************************************************
c     this subroutine has been modified to store all values of the
c     solution vector zy in every position eta (stored in zeig(j)).
c     the orthonormalization operation is done, at every step,
c     on all zeig(j) from j=1 to the present eta location.
c     Original coding by Anupa Bhajua. Modification by Marcio Mendonca
c************************************************************************

      amag1   = 0.d0
      amag2   = 0.d0
      amag3   = 0.d0
      gamma12 = dcmplx(0.d0)
      gamma13 = dcmplx(0.d0)
      gamma23 = dcmplx(0.d0)

c.......................find magnitude of the third vector (viscous soln)
      do i = 16, 21
        amag1 = amag1 + zy(i-12) * dconjg(zy(i-12))  
        amag2 = amag2 + zy(i-6)  * dconjg(zy(i-6))  
        amag3 = amag3 + zy(i)    * dconjg(zy(i))  
      end do
      amag1 = dsqrt(amag1)
      amag2 = dsqrt(amag2)
      amag3 = dsqrt(amag3)

      if ((max(amag1, amag2, amag3) .lt. 5.d1) .and.
     &    (jj .lt. jmax-1)) return

      amag1 = 0.d0
      amag2 = 0.d0

c..............................................normalise the third vector
      do i = 16, 21
        zy(i) = zy(i) / dcmplx(amag3)
      end do

c........................................................calculate gammas
      do i = 4, 9
        gamma13 = gamma13 + zy(i)   * dconjg(zy(i+12))       
        gamma23 = gamma23 + zy(i+6) * dconjg(zy(i+12))
      end do

c...........................orthogonalise the vectors and find magnitudes
c.....................................................second vector first
      do i = 10, 15
        zy(i) = zy(i) - gamma23 * zy(i+6)
        amag2 = amag2 + zy(i) * dconjg(zy(i))
      end do
      amag2 = dsqrt(amag2)

c.................................................normalize second vector
      do i = 10, 15
        zy(i) = zy(i) / dcmplx(amag2)
      end do

c..........................................orthogonalise the first vector
c...............................................use modified Gram-Schmidt 
c........................................a^1(1) = a(1) - (q(3)*a(1))*q(3)
c..............................a^2(1) = a(1) = a^1(1) - (q(2)*a^(1))*q(2)
      do i = 4, 9
        zy(i) =  zy(i) - gamma13 * zy(i+12)
      end do

      do i = 4, 9
        gamma12 = gamma12 + zy(i) * dconjg(zy(i+6))
      end do

c.......................orthogonalise the first vector and find magnitude
      do i = 4, 9
        zy(i) = zy(i) - gamma12 * zy(i+6)
        amag1 = amag1 + zy(i) * dconjg(zy(i))
      end do
      amag1 = dsqrt(amag1)

c..................................................normalise first vector
      do i = 4, 9
        zy(i) = zy(i) / dcmplx(amag1)
      end do

c.........................operate on all vector below the present vector
      do j = 1, jj-1
        do i = 16, 21
          zeig(i,j) = zeig(i,j) / dcmplx(amag3)   
        end do

        do i = 10, 15
          zeig(i,j) = zeig(i,j) - gamma23 * zeig(i+6,j)
          zeig(i,j) = zeig(i,j) / dcmplx(amag2)
        end do

        do i = 4, 9
          zeig(i,j) = zeig(i,j) - gamma13 * zeig(i+12,j)
     &                          - gamma12 * zeig(i+6,j)
          zeig(i,j) = zeig(i,j) / dcmplx(amag1)
        end do
      end do

      return
      end

**********************************************************************      

      subroutine out(eta, y, zfi, fo, h, hp, jedge)

      implicit none
      integer jedge, j, jurmax, juimax, jvrmax, jvimax, jwrmax, jwimax
      include 'par.f'

      complex*16 zu(4,jm), zfi(6,jm), alfa, alfa0, re, freq,
     &           ubar, ubarpp, z1, z2, z3
      real*8 u(jm), up(jm), upp(jm), ul(jm),
     &       fo(3,jm), eta(jm), y(jm), h(jm), hp(jm),
     &       f, fp, fpp, f3p, beta, re0, w0, m, beta_fs,
     &       urmax, uimax, vrmax, vimax, wrmax, wimax
      common /data/ re, alfa, alfa0, freq, beta, re0, w0, beta_fs

      urmax  = 0.d0
      uimax  = 0.d0
      vrmax  = 0.d0
      vimax  = 0.d0
      wrmax  = 0.d0
      wimax  = 0.d0
      jurmax = 0
      juimax = 0
      jvrmax = 0
      jvimax = 0
      jwrmax = 0
      jwimax = 0

      m = beta_fs / (2.d0 - beta_fs) 

      do j = 1, jedge-1
c........................................................stream function
        f   = fo(1,j)
        fp  = fo(2,j)
        fpp = fo(3,j)
        f3p = - ( m + 1.d0 ) / 2.d0 * fo(1,j) * fo(3,j)
     &        - m * ( 1.d0 - fo(2,j)**2 )    

c........................................mean flow velocity distribution
        u(j)   = fp
        up(j)  = fpp
        upp(j) = f3p

c...............................................................v = vbar
c..............................from continuity...ubar = vbar' * i / alfa
c.....................................alfa * ubar = alfa0 * u + beta * w
c.................................vorticity omega = beta * u - alfa0 * w
c.....................u = (beta * omega + alfa * alfa0 * ubar) / alfa**2
c....................................w = (alfa * ubar - alfa * u) / beta
c..............................................................zu(1) = u
c..............................................................zu(2) = v
c..............................................................zu(3) = w
c..............................................................zu(4) = p
        ubar    = zfi(2,j) * im / alfa
        ubarpp  = zfi(4,j) * im / alfa
        zu(1,j) = (dcmplx(beta) * zfi(5,j) + alfa * alfa0 * ubar)
     &          / alfa**2
        zu(2,j) = zfi(1,j)
        if (beta .eq. 0.d0) then
          zu(3,j) = dcmplx(0.d0)
         else
          zu(3,j) = (alfa * ubar - alfa0 * zu(1,j)) / dcmplx(beta)
        end if

c.....................................for pressure use equation for ubar
        z1      = im * alfa**2 / dcmplx(re) + freq - dcmplx(u(j)) * alfa
        z2      = im * dcmplx(up(j))
        z3      = im / re 
        zu(4,j) = (z1 * ubar + z2 * zfi(1,j) - z3 * ubarpp)
     &          * alfa0 / alfa**2
      end do

      do j = 1, jedge
        ul(j) = dreal(zu(1,j) * zexp(im * alfa0 * re0))
      end do

      zu(1,jedge) = dcmplx(0.d0)
      zu(2,jedge) = dcmplx(0.d0)
      zu(3,jedge) = dcmplx(0.d0)
      zu(4,jedge) = dcmplx(0.d0)

      open(unit = 11, file = 'uts.dat', status = 'unknown')
      open(unit = 12, file = 'vts.dat', status = 'unknown')
      open(unit = 13, file = 'wts.dat', status = 'unknown')
      open(unit = 14, file = 'pts.dat', status = 'unknown')

      do j = 1, jedge

c       write(11,*) j, zu(1,j)
c       write(12,*) j, zu(2,j)
c       write(13,*) j, zu(3,j)
c       write(14,*) j, zu(4,j)

c       if (dabs(dreal(zu(1,j))) .ge. dabs(urmax)) then
c         urmax  = dreal(zu(1,j))
c         jurmax = j
c       end if

c       if (dabs(dimag(zu(1,j))) .ge. dabs(uimax)) then
c         uimax  = dimag(zu(1,j))
c         juimax = j
c       end if

c       if (dabs(dreal(zu(2,j))) .ge. dabs(vrmax)) then
c         vrmax  = dreal(zu(2,j))
c         jvrmax = j
c       end if

c       if (dabs(dimag(zu(2,j))) .ge. dabs(vimax)) then
c         vimax  = dimag(zu(2,j))
c         jvimax = j
c       end if

c       if (dabs(dreal(zu(3,j))) .ge. dabs(wrmax)) then
c         wrmax  = dreal(zu(3,j))
c         jwrmax = j
c       end if

c       if ((dabs(dimag(zu(3,j))) .ge. dabs(wimax)) .and. (eta(j) .gt. 1.)) then
c         wimax  = dimag(zu(3,j))
c         jwimax = j
c       end if

      end do
      close(unit = 11)
      close(unit = 12)
      close(unit = 13)
      close(unit = 14)

c     do j = 2, jedge-1
c       write(10,110) zu(1,j), eta(j)
c       write(20,110) zu(2,j), eta(j)
c       write(30,110) zu(3,j), eta(j)
c       write(40,110) zu(4,j), eta(j)
c     end do
 110  format(3(e22.12))

c     open(unit = 15, file = 'alphats.dat', status = 'unknown')
c       write(15,*) im*alfa0
c     close(unit = 15)

c     open(unit = 10, file = 'ur.dat', access = 'append', status = 'unknown')
c     open(unit = 20, file = 'ui.dat', access = 'append', status = 'unknown')
c     open(unit = 30, file = 'vr.dat', access = 'append', status = 'unknown')
c     open(unit = 40, file = 'vi.dat', access = 'append', status = 'unknown')
c     open(unit = 50, file = 'wr.dat', access = 'append', status = 'unknown')
c     open(unit = 60, file = 'wi.dat', access = 'append', status = 'unknown')
c       write(10,*) w0/re0*1.d6, eta(jurmax)
c       write(20,*) w0/re0*1.d6, eta(juimax)
c       write(30,*) w0/re0*1.d6, eta(jvrmax)
c       write(40,*) w0/re0*1.d6, eta(jvimax)
c       write(50,*) w0/re0*1.d6, eta(jwrmax)
c       write(60,*) w0/re0*1.d6, eta(jwimax)
c       write(*,*) w0/re0*1.d6, jurmax, eta(jwrmax), re0
c     close(unit = 10)
c     close(unit = 20)
c     close(unit = 30)
c     close(unit = 40)
c     close(unit = 50)
c     close(unit = 60)

      return
      end

**********************************************************************      

      subroutine start(zy)

      implicit none
      include 'par.f'
      complex*16 zy(nm)
      real*8 wslope, etamax
      common/blas/ wslope, etamax

      zy(1)  = dcmplx(0.d0)
      zy(2)  = dcmplx(0.d0)
      zy(3)  = dcmplx(wslope)

      zy(4)  = dcmplx(0.d0)
      zy(5)  = dcmplx(0.d0)
      zy(6)  = dcmplx(1.d0)
      zy(7)  = dcmplx(0.d0)
      zy(8)  = dcmplx(0.d0)
      zy(9)  = dcmplx(0.d0)

      zy(10) = dcmplx(0.d0)
      zy(11) = dcmplx(0.d0)
      zy(12) = dcmplx(0.d0)
      zy(13) = dcmplx(1.d0)
      zy(14) = dcmplx(0.d0)
      zy(15) = dcmplx(0.d0)

      zy(16) = dcmplx(0.d0)
      zy(17) = dcmplx(0.d0)
      zy(18) = dcmplx(0.d0)
      zy(19) = dcmplx(0.d0)
      zy(20) = dcmplx(0.d0)
      zy(21) = dcmplx(1.d0)

      return
      end

**********************************************************************      

      subroutine rkdumb(zystart, nvar, zy, zeig, jmax, eta)

      implicit none
      integer j, jmax, i, nvar
      include 'par.f'
      complex*16 zystart(nvar), zy(nvar), zdydx(nm),
     &           zeig(nvar,jm), zyout(nm)
      real*8 eta(jm), h, x

      do i = 1,nvar
        zy(i)     = zystart(i)
        zeig(i,1) = zy(i)
      end do

      do j = 1, jmax-1
        call derivs(x, zy, zdydx)
        h = eta(j+1) - eta(j)
        call rk4(zy, zdydx, nvar, h, zyout)

        do i = 1, nvar
          zy(i) = zyout(i)
        end do

        call ortho(zy, zeig, j+1, jmax)

        do i = 1, nvar
          zeig(i,j+1) = zy(i)
        end do
      end do

      return
      end 

**********************************************************************      

      subroutine rk4(zy, zdydx, n, h, zyout)

      implicit none
      integer i, n
      include 'par.f'
      complex*16 zy(n), zdydx(n), zyout(n), zyt(nm),
     &           zdyt(nm),zdym(nm)
      real*8 h, hh, h6, x

      hh = h * 0.5d0
      h6 = h / 6.d0

      do i = 1, n
        zyt(i) = zy(i) + dcmplx(hh) * zdydx(i)
      end do

      call derivs(x, zyt, zdyt) 
      do i = 1, n
        zyt(i) = zy(i) + dcmplx(hh) * zdyt(i)
      end do

      call derivs(x, zyt, zdym)
      do i = 1, n
        zyt(i)  = zy(i)   + dcmplx(h) * zdym(i)
        zdym(i) = zdyt(i) + zdym(i)
      end do

      call derivs(x, zyt, zdyt)
      do i = 1, n
        zyout(i) = zy(i) + dcmplx(h6)
     &           * (zdydx(i) + zdyt(i) + dcmplx(2.d0) * zdym(i))
      end do

      return
      end
