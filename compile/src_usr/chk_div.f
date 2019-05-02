      subroutine chk_div(vx,vy,vz,divv)
C--------------------------------------------------------------------
C     Author: Prabal Negi
C     Check the divergence of the field
C
C     This routine is a modified version of the nek routine: chkptol
C 
C--------------------------------------------------------------------

      implicit none

      include 'SIZE_DEF'
      include 'SIZE'
      include 'MASS_DEF'
      include 'MASS'
      include 'INPUT_DEF'
      include 'INPUT'
      include 'TSTEP_DEF'
      include 'TSTEP'


      real VX(LX1,LY1,LZ1,LELT)
      real VY(LX1,LY1,LZ1,LELT)
      real VZ(LX1,LY1,LZ1,LELT)

      real DIVV (LX2,LY2,LZ2,LELT)                    ! divergence without mass
      real BDIVV(LX2,LY2,LZ2,LELT)                    ! Divergence with mass
      COMMON /SCRUZ/ BDIVV

      integer NTOT2
      real DNORM
      real DMAX

      real GLAMAX             ! function. Global absolute max
      real GLSC2              ! L2 inner product

      NTOT2 = NX2*NY2*NZ2*NELV
      CALL OPDIV (BDIVV,VX,VY,VZ)                     ! Has factor BM2
      CALL COL3 (DIVV,BDIVV,BM2INV,NTOT2)             ! Remove BM2 factor

      DNORM = SQRT(GLSC2(DIVV,BDIVV,NTOT2)/VOLVM2)    ! norm
      DMAX = GLAMAX(DIVV,NTOT2)                       ! absolute max
C
      if (nio.eq.0) WRITE (6,101) istep, time, 'DNORM, DMAX',DNORM,DMAX
  101 format((I6,2x,E15.8E3,2x,A11,2x,2(E15.8E3,2x)))   
C
      RETURN
      END

