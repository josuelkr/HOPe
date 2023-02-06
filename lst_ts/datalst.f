      program geradatalst

      implicit none
      integer f
      include '../par.for'

      write(*,*) 'This analysis is for a 2D (0) or a 3D (1) T-S wave ?'
      read(*,*) f
      open(2, file = 'datalst.dat', form = 'formatted')
        write(2,*) dsqrt(Re), omega/dsqrt(Re), beta*dble(f)/dsqrt(Re),
     &             imax, dx
      close(unit = 2)

      stop
      end
