!====================================================================== 
!     Author: Prabal Negi
!     Description: Wrapper routines for Lapack
!     List of Wrapper Routines:
!     1) wrp_lls:       dgels: Linear Least Squares
!     2) wrp_svd:       dgesdd: SVD decomposition 
!     3) wrp_rschur:    dgees: Real Schur Decomposition
!     4) wrp_cschur:    cgees: Complex (Single) Schur Decomposition
!     5) wrp_zschur:    zgees: Complex (Double) Schur Decomposition           
!     6) wrp_zgeinv:    zgetrf->zgetri: Complex (Double) Matrix Inversion
!     7) wrp_zgeeig:    zgeev: Complex (Double) Matrix Eigenvalue Decomposition
!     8) wrp_dgeinv:    dgetrf->dgetri: Real (Double) Matrix Inversion
!     9) wrp_dgeev:     dgeev: Real (Double) Matrix Eigenvalue Decomposition      
!======================================================================       
!---------------------------------------------------------------------- 
      subroutine wrp_lls(m,n,nrhs,Amat,lda,rhs,ldb)

!     LAPACK interface for Linear Least Squares.

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'
      
      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  DBDSDC did not converge, updating process failed.                          
      character*1 trans ! Specifies which problem to solve:
                        ! = 'N': LS solution of the problem Ax-b=0;
                        ! = 'T': LS solution of the problem (A^T)x-b=0;
                        ! = 'C': LS solution of the problem (A^H)x-b=0;

      integer m         ! no rows in matrix
      integer n         ! no columns in matrix
      integer lda       ! leading dimension of the matrix
      integer ldb       ! leading dimension of rhs vector/matrix
      integer nrhs      ! no of RHS vectors

      real Amat(lda,n)        ! A(lda,n)
      real rhs(ldb,nrhs)      ! b(ldb,nrhs)


!     Define LAPACK-variables
      trans = 'N'

!     Compute Least-Squares in double precision
      call dgels(trans,m,n,nrhs,Amat,lda,rhs,ldb,LLS_WKR,LLS_WRL
     $           ,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: the i:th argment had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: DBDSDC did not converge, updating process failed.'
     $        , info
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGELS: successful exit!'
         if (nid.eq.0) write(6,*) '         Optimal LWKR=',
     $        int(LLS_WKR(1)), LLS_WRL
      endif
            
      return
      end subroutine wrp_lls
!----------------------------------------------------------------------
      subroutine wrp_svd(m,n,Amat,lda,Sigma,U,ldu,VT,ldvt)

!     LAPACK interface for SVD.

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'
      
      character*1 jobz  !  Specifies options for computing all or part of the matrix U:
                        ! 'A':  all M columns of U and all N rows of V**T are
                        !    returned in the arrays U and VT;
                        ! 'S':  the first min(M,N) columns of U and the first
                        !   min(M,N) rows of V**T are returned in the arrays U and VT;
                        ! 'O':  If M >= N, the first N columns of U are overwritten
                        !   on the array A and all rows of V**T are returned in the array VT;
                        !   otherwise, all columns of U are returned in the array U 
                        !   and the first M rows of V**T are overwritten in the array A;
                        ! 'N':  no columns of U or rows of V**T are computed.

      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  DBDSDC did not converge, updating process failed.                          
      
      integer m         ! no rows in matrix
      integer n         ! no columns in matrix
      integer lda       ! leading dimension of the matrix
      integer ldu       ! leading dimension of u
      integer ldvt      ! leading dimension of VT

      real U(ldu,n)     ! Left Singular vectors       ! not referenced
      real Amat(lda,n)  ! compute SVD of Amat
      real Sigma(n)     ! Singular values
      real VT(ldvt,n)   ! Right singular vectors

      integer i

      character*23 fname

!     Define LAPACK-variables
      jobz = 'O'                ! return the n first singular vectors 

!     Compute SVD in double precision with divide-and-conquer
      call dgesdd(jobz,m,n,Amat,lda,Sigma,U,ldu,VT,ldvt,
     $     SVD_WKR,SVD_WRL,SVD_WKI,info)

!     Error-check
      if (info.lt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: the i:th argment had an illegal value.', abs(info)
         call exitt
      elseif (info.gt.0) then
         if (nid.eq.0) write(6,*)
     $        'ERROR: DBDSDC did not converge, updating process failed.'
     $        , info
         call exitt
      else
         if (nid.eq.0) write(6,*) 'DGESDD: successful exit!'
         if (nid.eq.0) write(6,*) '         Optimal LWKR=',
     $        int(SVD_WKR(1)), SVD_WRL
      endif
            
!     Output singular values and estimated error bounds
      write(fname,'(A19)') 'singular_values.txt'
      if (nid.eq.0) then
         open(61,file=fname,action='write',status='unknown')
         write(61,'(TR2,A2,TR18,A6)')
     $        'I;','sigma;'

         do i=1,min(m,n)
            write(61,'(I4,4G24.16)') i, SIGMA(i)
         enddo
         close(61)
      endif
      
      return
      end subroutine wrp_svd
!----------------------------------------------------------------------
      subroutine wrp_rschur(Amat,lda,n,VS,ldvs,wr,wi)

!     LAPACK interface for Real Schur decomposition.

      implicit none

      include 'SIZE'
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

      real Amat(lda,n)  ! compute Schur form of Amat
      real wr(n)        ! diagonal values? (real)
      real wi(n)        ! diagonal values? (imaginary)
      real VS(ldvs,n)   ! Orthnormal Schur Vectors

      integer i

!     Define LAPACK-variables
      jobvs = 'V'
      esort = 'N'

!     Compute real-Schur decomposition in double precision
      call dgees(jobvs,esort,selec,n,Amat,lda,sdim,wr,wi,VS,ldvs,
     $     RSCHUR_WKR,RSCHUR_WRL,RSCHUR_WKB,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
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
        endif               
        call exitt
      else
         if (nid.eq.0) write(6,*) 'DGEES: successful exit!'
         if (nid.eq.0) write(6,*) '         Optimal LWKR=',
     $        int(RSCHUR_WKR(1)), RSCHUR_WRL
      endif
            
      return
      end subroutine wrp_rschur
!----------------------------------------------------------------------
      subroutine wrp_cschur(Amat,lda,n,VS,ldvs,w)

!     LAPACK interface for Complex Schur decomposition.

      implicit none

      include 'SIZE'
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

!     Define LAPACK-variables
      jobvs = 'N'
      esort = 'N'

!     Compute complex-Schur decomposition in Single precision
      call cgees(jobvs,esort,selec,n,Amat,lda,sdim,w,VS,ldvs,
     $     CSCHUR_WKC,CSCHUR_WCL,CSCHUR_WKR,CSCHUR_WKB,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
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
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'CGEES: successful exit!'
        if (nid.eq.0) write(6,*) '         Optimal WCL=',
     $       int(CSCHUR_WKC(1)), CSCHUR_WCL
      endif
            
      return
      end subroutine wrp_cschur
!----------------------------------------------------------------------
      subroutine wrp_zschur(Amat,lda,n,VS,ldvs,w)

!     LAPACK interface for Complex Schur decomposition.

      implicit none

      include 'SIZE'
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

!     Define LAPACK-variables
      jobvs = 'V'
      esort = 'N'


!     Compute complex-Schur decomposition in double precision
      call zgees(jobvs,esort,selec,n,Amat,lda,sdim,w,VS,ldvs,
     $     ZSCHUR_WKC,ZSCHUR_WCL,ZSCHUR_WKR,ZSCHUR_WKB,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
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
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'ZGEES: successful exit!'
        if (nid.eq.0) write(6,*) '         Optimal WCL=',
     $       int(ZSCHUR_WKC(1)), ZSCHUR_WCL
      endif
            
      return
      end subroutine wrp_zschur
!----------------------------------------------------------------------

      subroutine wrp_zgeinv(Amat,lda,n)

!     LAPACK interface for general Complex matrix inversion.

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'

      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  See Doc for error 

      integer m         ! rows of the Matrix A(lda,n)
      integer n         ! Columns of the Matrix A(lda,n)
      integer lda       ! leading dimension of the matrix A

      complex Amat(lda,n)  ! compute Schur form of Amat


!     Perform LU decomposition
      m=n
      call zgetrf(m,n,Amat,lda,ZGEINV_WKI,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: U(i,i) is exactly zero in zgetrf.'
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'ZGETRF: successful exit!'
      endif

!     Matrix inversion
      call zgetri(n,Amat,lda,ZGEINV_WKI,ZGEINV_WKC,ZGEINV_WCL,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: U(i,i) is exactly zero in zgetri.'
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'ZGETRI: successful exit!'
      endif


      return
      end subroutine wrp_zgeinv
!----------------------------------------------------------------------       

      subroutine wrp_zgeeig(Amat,lda,n,w,VL,ldvl,VR,ldvr)

!     LAPACK interface for Eigenpair decomposition of a complex matrix.

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'

      character jobvl*1 ! = 'N': Left eigenvectors are not computed
                        ! = 'V': Compute Left eigenvectors

      character jobvr*1 ! = 'N': Right eigenvectors are not computed
                        ! = 'V': Compute Right eigenvectors


      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  See Doc for error. 

      integer n         ! Order of the Matrix A(lda,n)
      integer lda       ! leading dimension of the matrix A
      integer ldvl      ! leading dimension of the matrix VL(ldvl,n)
      integer ldvr      ! leading dimension of the matrix VR(ldvr,n)

      complex Amat(lda,n)  ! compute eigenpairs of Amat
      complex VL(ldvl,n)  ! Left Eigenvectors
      complex VR(ldvr,n)  ! Right Eigenvectors
     
      complex w(n)      ! eigenvalues

      

!     Perform Eigen decomposition
      jobvl = 'N'
      jobvr = 'N'
      call zgeev(jobvl,jobvr,n,Amat,lda,w,VL,ldvl,VR,ldvr,ZEIG_WKC,
     $           ZEIG_WCL,ZEIG_WKR,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: QR algorithm failed to compute all eigenvalues. ',
     $     'No eigenvectors computed.'          
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'ZGEEV: successful exit!'
        if (nid.eq.0) write(6,*) '           Optimal WCL=',
     $       int(ZEIG_WKC(1)), ZEIG_WCL
      endif

      return
      end subroutine wrp_zgeeig
!----------------------------------------------------------------------       

      subroutine wrp_dgeinv(Amat,lda,n)

!     LAPACK interface for general Real (Double) matrix inversion.

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'

      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  See Doc for error 

      integer m         ! rows of the Matrix A(lda,n)
      integer n         ! Columns of the Matrix A(lda,n)
      integer lda       ! leading dimension of the matrix A

      real Amat(lda,n)  ! compute Schur form of Amat


!     Perform LU decomposition
      m=n
      call dgetrf(m,n,Amat,lda,DGEINV_WKI,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: U(i,i) is exactly zero in dgetrf.'
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'DGETRF: successful exit!'
      endif

!     Matrix inversion
      call dgetri(n,Amat,lda,DGEINV_WKI,DGEINV_WKR,DGEINV_WRL,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: U(i,i) is exactly zero in dgetri.'
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'DGETRI: successful exit!'
      endif


      return
      end subroutine wrp_dgeinv
!----------------------------------------------------------------------       

      subroutine wrp_dgeev(Amat,lda,n,wr,wi,VL,ldvl,VR,ldvr)

!     LAPACK interface for Eigenpair decomposition of a general matrix.
!     Double precision        

      implicit none

      include 'SIZE'
      include 'WRP_LAPACK'

      character jobvl*1 ! = 'N': Left eigenvectors are not computed
                        ! = 'V': Compute Left eigenvectors

      character jobvr*1 ! = 'N': Right eigenvectors are not computed
                        ! = 'V': Compute Right eigenvectors


      integer info      ! = 0:  successful exit.
                        ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                        ! > 0:  if INFO = i, the QR algorithm failed to compute all the
                        !       eigenvalues, and no eigenvectors have been computed;
                        !       elements i+1:N of WR and WI contain eigenvalues which
                        !       have converged.

      integer n         ! Order of the Matrix A(lda,n)
      integer lda       ! leading dimension of the matrix A
      integer ldvl      ! leading dimension of the matrix VL(ldvl,n)
      integer ldvr      ! leading dimension of the matrix VR(ldvr,n)

      complex Amat(lda,n)  ! compute eigenpairs of Amat
      complex VL(ldvl,n)  ! Left Eigenvectors
      complex VR(ldvr,n)  ! Right Eigenvectors
     
      real wr(n)          ! eigenvalues real part
      real wi(n)          ! eigenvalues imaginary part

      if (SW_DGEEV.eq.-1) then
        if (nio.eq.0) write(6,*) 
     $    'WRP_DGEEV not switched on in WRP_LAPACK. SW_DGEEV=',SW_DGEEV
        call exitt
      endif  

!     Perform Eigen decomposition
      jobvl = 'N'
      jobvr = 'V'
      call dgeev(jobvl,jobvr,n,Amat,lda,wr,wi,VL,ldvl,VR,ldvr,
     $           DGEEV_WKR,DGEEV_WRL,info)

!     Error-check
      if (info.lt.0) then
        if (nid.eq.0) write(6,*)
     $       'ERROR: the i:th argment had an illegal value.', abs(info)
        call exitt
      elseif (info.gt.0) then
        if (nid.eq.0) then
          write(6,*)
     $     'ERROR: QR algorithm failed to compute all eigenvalues. ',
     $     'No eigenvectors computed.'          
     $     ,info
        endif               
        call exitt
      else
        if (nid.eq.0) write(6,*) 'DGEEV: successful exit!'
        if (nid.eq.0) write(6,*) 'Optimal WRL=',
     $       int(DGEEV_WKR(1)), DGEEV_WRL
      endif

      return
      end subroutine wrp_dgeev
!----------------------------------------------------------------------       








