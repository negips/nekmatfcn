!====================================================================== 
!     Author: Prabal Negi
!     Description: Wrapper routines for Lapack
!
!======================================================================       
!----------------------------------------------------------------------
      subroutine wrp_cschur(Amat,lda,n,VS,ldvs,w)

!     LAPACK interface for Complex Schur decomposition.

      implicit none

!      include 'SIZE_DEF'
!      include 'SIZE'
      include 'WRP_LAPACK'
      
      character*1 jobvs ! 'N': Schur vectors are not computed
                        ! 'V': Schur vectors are computed

      character*1 esort ! 'N': unsorted Eigenvalues
                        ! 'S': Ordered Eigenvalues.
                        ! Currently hard-coded to 'N'.
                        ! 'S' needs a reverse communication routine.

      logical selec     ! Not used                        


      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  DBDSDC did not converge, updating process failed.                          
      integer sdim      ! not used if esort='N'

      integer n         ! order of the Matrix A(n,n)
      integer lda       ! leading dimension of the matrix A
      integer ldvs      ! leading dimension of VS

      complex Amat(lda,n)  ! compute Schur form of Amat
      complex w(n)         ! Eigenvalues
      complex VS(ldvs,n)   ! Orthnormal Schur Vectors

      integer i

      character*23 fname

!     Define LAPACK-variables
      jobvs = 'N'
      esort = 'N'


!     Compute complex-Schur decomposition in double precision
      call cgees(jobvs,esort,selec,n,Amat,lda,sdim,w,VS,ldvs,
     $     CSCHUR_WORKC,CSCHUR_WCL,CSCHUR_WORKR,CSCHUR_WORKB,info)

!     Error-check
      if (info.lt.0) then
        write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
!        call exitt
      elseif (info.gt.0) then
          if (info.le.n) then
            write(6,*)
     $       'ERROR: QR Algorithm failed to compute all eigenvalues'
     $       ,info
          elseif (info.eq.(n+1)) then
            write(6,*)
     $       'ERROR: Eigenvalues could not be ordered.'
     $       ,info
          elseif (info.eq.(n+2)) then
            write(6,*)
     $       'ERROR: Some round-off nonesense. See Doc. pg 222'
     $       ,info
          endif
!        call exitt
      else
            write(6,*) 'CGEES: successful exit!'
            write(6,*) '         Optimal WCL=',
     $       int(CSCHUR_WORKC(1)), CSCHUR_WCL
      endif
            
      return
      end subroutine wrp_cschur
!----------------------------------------------------------------------



